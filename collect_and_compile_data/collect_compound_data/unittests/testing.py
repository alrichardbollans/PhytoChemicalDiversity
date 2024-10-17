import pandas as pd
import unittest
from collect_and_compile_data.collect_compound_data import NP_PATHWAYS
from collect_and_compile_data.collect_compound_data import taxon_group_resolution_methods


class TestTaxonGroupResolutionMethods(unittest.TestCase):

    def setUp(self):
        self.data = {
            'Genus': ['Terp_fat', 'None', 'Terp_fat_poly', 'Terp_poly'],
            'Terpenoids': [1, 0, 1, 1],
            'Fatty_acids': [1, 0, 1, 0],
            'Polyketides': [0, 0, 1, 1],
            'Carbohydrates': [0, 0, 0, 0],
            'Amino_acids_and_Peptides': [0, 0, 0, 0],
            'Shikimates_and_Phenylpropanoids': [0, 0, 0, 0],
            'Alkaloids': [0, 0, 0, 0],
        }

    def test_split_multiple_pathways_into_duplicate_rows(self):
        df = pd.DataFrame(self.data)
        result = taxon_group_resolution_methods.split_multiple_pathways_into_duplicate_rows(df)

        # Assert that the resulting dataframe has more rows than the original
        self.assertGreater(len(result), len(df))

        # Assert that each row in the resulting dataframe has a sum <= 1
        for index, row in result.iterrows():
            self.assertLessEqual(row[NP_PATHWAYS].sum(), 1)


if __name__ == '__main__':
    unittest.main()
