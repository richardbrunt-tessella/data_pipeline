import unittest

from mrtarget.modules.Association import AssociationMap


class AssociationMapTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.as_map = AssociationMap(id='testid')
        cls.as_map.add_datasource('chembl', 'test1')
        cls.as_map.add_datasource('chembl', 'test2')
        cls.as_map.add_datasource('chembl', 'test3')
        cls.as_map.add_datasource('eva_somatic', 'test3')
        cls.as_map.add_datasource('eva_somatic', 'test1')
        cls.as_map.add_datasource('cancer_gene_census', 'test3')

    def test_encode_vector(self):
        encoded_vector = AssociationMap.encode_vector(['1', '3', '5'], map(str, range(10)))
        self.assertEquals(encoded_vector, ' '.join(['n0', '1', 'n2', '3', 'n4', '5', 'n6', 'n7', 'n8', 'n9']))

    def test_map_creation(self):
        as_map = self.as_map
        self.assertEquals(len(as_map.overall), 3)
        self.assertEquals(len(as_map.datatype['known_drug']), 3)
        self.assertEquals(len(as_map.datatype['somatic_mutation']), 2)
        self.assertEquals(len(as_map.datasource['chembl']), 3)
        self.assertEquals(len(as_map.datasource['eva_somatic']), 2)
        self.assertEquals(len(as_map.datasource['cancer_gene_census']), 1)
        json_as_map = as_map.to_json()
        self.assertTrue(isinstance(json_as_map, str))
        self.assertIn('test1', json_as_map)
        self.assertIn('test2', json_as_map)
        self.assertIn('test3', json_as_map)
        self.assertIn('testid', json_as_map)

    def test_encode_all_vectors(self):
        self.as_map.encode_all_vectors(['test0',
                                        'test1',
                                        'test2',
                                        'test3',
                                        'test4',
                                        'test5', ])
        self.assertEquals(self.as_map.overall_vector, ' '.join(['ntest0',
                                                         'test1',
                                                         'test2',
                                                         'test3',
                                                         'ntest4',
                                                         'ntest5', ]))
