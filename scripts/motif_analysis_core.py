"""MEME/TOMTOM summaries and independent AME analysis for archived CHO promoters."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal
import json
import hashlib
from pathlib import Path
import subprocess
import tarfile
import tempfile
import xml.etree.ElementTree as ET
import pandas as pd

DATABASES = {'TF_AllVertebrate': 'tomtom_all_out', 'TF_Mouse': 'tomtom_mmus_out'}
COLUMNS = ['iModulon_ID','motif_name','motif_id','Total_Sites','TF','TOMTOM_qvalue','MEME_Evalue','TF_Confidence','Source']
AME_COLUMNS = ['iModulon_ID','motif_ID','motif_alt_ID','consensus','AME_pvalue','AME_adj_pvalue','AME_Evalue','TP','TP_percent','FP','FP_percent','pos','neg','AME_Significance','AME_Evalue_significant']

def table(path):
    if not path.exists():
        return None
    try:
        return pd.read_csv(path, sep='\t', comment='#', dtype=str)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

def unpack(archive, destination):
    manifest = archive.with_name(archive.name.replace('_inputs.tar.gz', '_input_manifest.json'))
    if manifest.exists():
        expected = json.loads(manifest.read_text())['sha256']
        with archive.open('rb') as stream:
            actual = hashlib.file_digest(stream, 'sha256').hexdigest()
        if actual != expected:
            raise ValueError(f'Archive checksum mismatch: {archive}')
    with tarfile.open(archive) as tar:
        tar.extractall(destination, filter='data')

def summarize(root, output, prefix, excluded_modules=None):
    excluded_modules = excluded_modules or {}
    output.mkdir(parents=True, exist_ok=True)
    rows, ame_rows, coverage, inventory = [], [], [], []
    for folder in sorted((p for p in root.iterdir() if p.is_dir() and p.name.isdigit()), key=lambda p:int(p.name)):
        module = 'iModulon_' + folder.name
        xml = folder/'meme_out/meme.xml'
        status = {'iModulon_ID':module, 'MEME':'missing', 'TF_Mouse':'missing', 'TF_AllVertebrate':'missing', 'AME':'missing'}
        if int(folder.name) in excluded_modules:
            status.update({key: 'excluded_single_gene' for key in ['MEME','TF_Mouse','TF_AllVertebrate','AME']})
            coverage.append(status)
            continue
        motifs = {}
        if xml.exists() and xml.stat().st_size:
            for m in ET.parse(xml).findall('.//motifs/motif'):
                rec = {'iModulon_ID':module, 'motif_name':m.get('alt'), 'motif_id':m.get('name'), 'MEME_Evalue':m.get('e_value')}
                inventory.append(rec)
                if Decimal(m.get('e_value')) < Decimal('0.01'):
                    motifs[m.get('name')] = rec
            status['MEME'] = 'complete'
        counts = {}
        countfile = folder/'fimo_site_counts.csv'
        if countfile.exists():
            counts = pd.read_csv(countfile).set_index('motif_id')['Total_Sites'].to_dict()
        for source, subdir in DATABASES.items():
            hits = table(folder/subdir/'tomtom.tsv')
            if hits is None:
                continue
            status[source] = 'complete' if len(hits) else 'no_hits'
            labels = {}
            labelpath = folder/subdir/'tomtom.xml'
            if labelpath.exists():
                for target in ET.parse(labelpath).findall('.//targets/motif'):
                    ident = target.get('id')
                    labels[ident] = (target.get('alt') or ident).removeprefix(ident + '.')
            for hit in hits.to_dict('records'):
                if hit['Query_ID'] not in motifs:
                    continue
                ident = hit['Target_ID']
                q = float(hit['q-value'])
                rows.append(dict(motifs[hit['Query_ID']], Total_Sites=counts.get(hit['Query_ID']), TF=labels.get(ident, ident), TOMTOM_qvalue=q, TF_Confidence='HIGH' if q < .05 else 'LOW', Source=source, Target_ID=ident))
        ame = table(folder/'ame_out/ame.tsv')
        if ame is not None:
            status['AME'] = 'complete' if len(ame) else 'no_hits'
            mapping={'p-value':'AME_pvalue','adj_p-value':'AME_adj_pvalue','E-value':'AME_Evalue','%TP':'TP_percent','%FP':'FP_percent'}
            for hit in ame.rename(columns=mapping).to_dict('records'):
                hit.update(iModulon_ID=module, AME_Significance='Significant' if Decimal(hit['AME_adj_pvalue']) < Decimal('.05') else 'Not significant', AME_Evalue_significant=Decimal(hit['AME_Evalue']) < Decimal('.05'))
                ame_rows.append(hit)
        coverage.append(status)
    detailed = pd.DataFrame(rows, columns=COLUMNS+['Target_ID'])
    detailed = detailed.sort_values(['TOMTOM_qvalue','Target_ID'], kind='stable').drop_duplicates(['iModulon_ID','motif_id','Source','Target_ID'])
    summary = detailed.drop_duplicates(['iModulon_ID','motif_id','Source','TF'])[COLUMNS]
    summary.to_csv(output/f'{prefix}_meme_tomtom_summary.csv', index=False)
    detailed.to_csv(output/f'{prefix}_tomtom_detailed_hits.csv', index=False)
    for source in DATABASES:
        summary.loc[summary.Source == source].to_csv(output/f'{prefix}_{source}_summary.csv', index=False)
    pd.DataFrame(ame_rows, columns=AME_COLUMNS).to_csv(output/f'{prefix}_ame_summary.csv', index=False)
    pd.DataFrame(inventory, columns=['iModulon_ID','motif_name','motif_id','MEME_Evalue']).to_csv(output/f'{prefix}_meme_motif_inventory.csv', index=False)
    status = pd.DataFrame(coverage)
    status.to_csv(output/f'{prefix}_motif_coverage.csv', index=False)
    missing = status.loc[(status.iloc[:,1:] == 'missing').any(axis=1)]
    print(f'{prefix}: {len(status)} modules; {len(summary)} TF matches; {len(ame_rows)} AME rows; {len(missing)} modules with missing outputs.')
    return summary, status

def excluded_modules(location, cohort):
    manifest = json.loads((location/f'{cohort}_motif_input_manifest.json').read_text())
    return {record['iModulon']: record['analysis_exclusion'] for record in manifest['modules'] if record.get('analysis_exclusion')}

def summary_cli(cohort, location):
    parser = argparse.ArgumentParser(description='Summarize MEME XML and TOMTOM results; AME is read independently. No FIMO, CentriMo or SpaMo required.')
    parser.add_argument('--results-dir', type=Path)
    parser.add_argument('--output-dir', type=Path, default=location/'results')
    args = parser.parse_args()
    if args.results_dir:
        summarize(args.results_dir, args.output_dir, cohort, excluded_modules(location, cohort))
    else:
        with tempfile.TemporaryDirectory() as tmp:
            unpack(location/f'{cohort}_motif_inputs.tar.gz', tmp)
            summarize(Path(tmp), args.output_dir, cohort, excluded_modules(location, cohort))

def run_cli(cohort, location):
    parser = argparse.ArgumentParser(description='Rerun MEME, two TOMTOM databases and AME on the supplied original promoter sequences (MEME Suite required).')
    db_root = Path(__file__).resolve().parents[1]/'databases/JASPAR2024'
    parser.add_argument('--vertebrate-db', type=Path, default=db_root/'JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.meme')
    parser.add_argument('--mouse-db', type=Path, default=db_root/'JASPAR2024_CORE_non-redundant_pfms_meme_musmusculus.meme')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=1)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--modules', nargs='+', type=int)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    if args.jobs < 1 or args.threads < 1:
        parser.error('jobs and threads must be positive')
    for db in [args.vertebrate_db, args.mouse_db]:
        if not db.is_file():
            parser.error(f'Database not found: {db}')
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error('Use a new or empty output directory to preserve previous results.')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        unpack(location/f'{cohort}_motif_inputs.tar.gz', tmp)
        folders = sorted((p for p in Path(tmp).iterdir() if p.is_dir() and p.name.isdigit()), key=lambda p:int(p.name))
        if args.modules:
            available = {int(p.name) for p in folders}
            if set(args.modules) - available:
                parser.error('Requested modules are absent from the archive')
            folders = [p for p in folders if int(p.name) in args.modules]
        exclusions = excluded_modules(location, cohort)
        def worker(source):
            import shutil
            dest = args.output_dir.resolve()/source.name
            dest.mkdir()
            if int(source.name) in exclusions:
                (dest/'analysis_exclusion.txt').write_text(exclusions[int(source.name)] + '\n')
                return
            for name in ['promoters_primary.fa','promoters_control.fa']:
                shutil.copy2(source/name, dest/name)
            primary, control = str(dest/'promoters_primary.fa'), str(dest/'promoters_control.fa')
            commands = [
                ['fasta-get-markov','-m','1',primary,str(dest/'promoters.bg')],
                ['meme',primary,'-dna','-oc',str(dest/'meme_out'),'-nostatus','-mod','zoops','-nmotifs','10','-minw','6','-maxw','40','-revcomp','-bfile',str(dest/'promoters.bg'),'-p',str(args.threads)],
                ['tomtom','-oc',str(dest/'tomtom_all_out'),'-dist','pearson',str(dest/'meme_out/meme.txt'),str(args.vertebrate_db.resolve())],
                ['tomtom','-oc',str(dest/'tomtom_mmus_out'),'-dist','pearson',str(dest/'meme_out/meme.txt'),str(args.mouse_db.resolve())],
                ['ame','--oc',str(dest/'ame_out'),'--control',control,'--scoring','avg','--method','fisher',primary,str(args.vertebrate_db.resolve())]]
            (dest/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
            if not args.dry_run:
                with (dest/'execution.log').open('w') as log:
                    for command in commands:
                        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            list(pool.map(worker, folders))
    if not args.dry_run:
        summarize(args.output_dir, args.output_dir/'summary', cohort, excluded_modules(location, cohort))
    print('Command preparation complete.' if args.dry_run else 'Motif analysis complete.')
