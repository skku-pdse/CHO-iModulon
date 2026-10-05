import sys
from pathlib import Path
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.motif_analysis_core import summarize

class MotifSummaryTests(unittest.TestCase):
    def test_database_separation_without_downstream_tools(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); module=root/'0'; (module/'meme_out').mkdir(parents=True)
            (module/'meme_out/meme.xml').write_text('<MEME><motifs><motif name="ACGT" alt="MEME-1" e_value="1e-400" /></motifs></MEME>')
            for directory,ids in [('tomtom_all_out',['A','B','C']),('tomtom_mmus_out',['A','A'])]:
                p=module/directory;p.mkdir()
                (p/'tomtom.tsv').write_text('Query_ID\tTarget_ID\tq-value\n'+''.join(f'ACGT\t{x}\t0.01\n' for x in ids))
            (module/'ame_out').mkdir();(module/'ame_out/ame.tsv').write_text('# no significant hits\n')
            frame,coverage=summarize(root,root/'summary','test')
            self.assertEqual(len(frame),4)
            self.assertTrue(frame.Total_Sites.isna().all())
            self.assertEqual(coverage.iloc[0].AME,'no_hits')
            self.assertEqual(frame.iloc[0].MEME_Evalue,'1e-400')
    def test_missing_module_is_not_negative_result(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'106').mkdir()
            frame,coverage=summarize(root,root/'summary','test')
            self.assertTrue(frame.empty)
            self.assertEqual(coverage.iloc[0].MEME,'missing')
            self.assertEqual(coverage.iloc[0].AME,'missing')

    def test_single_gene_exclusion_is_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'106').mkdir()
            frame,coverage=summarize(root,root/'summary','expanded', {106:'Single-gene iModulon'})
            self.assertTrue(frame.empty)
            self.assertTrue((coverage.iloc[0,1:]=='excluded_single_gene').all())

if __name__=='__main__':unittest.main()
