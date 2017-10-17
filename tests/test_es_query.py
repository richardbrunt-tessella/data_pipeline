import unittest

from mrtarget.Settings import Config
from mrtarget.common.ElasticsearchQuery import ESQuery


class AssociationMapTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.es_query = ESQuery()

    @unittest.skipIf((not Config.RELEASE_VERSION) or (not Config.ELASTICSEARCH_NODES),
                     "not supported in this library version")
    def test_distinct_targets(self):
        total_targets = self.es_query.count_all_targets()
        targets_with_evidence = list(self.es_query.get_all_target_ids_with_evidence_data())
        self.assertLess(len(targets_with_evidence), total_targets)
        self.assertTrue(targets_with_evidence[0].startswith('ENSG'))

    def test_distinct_disease(self):
        total_diseases = self.es_query.count_all_diseases()
        diseases_with_evidence = list(self.es_query.get_all_disease_ids_with_evidence_data())
        self.assertLess(len(diseases_with_evidence), total_diseases)
