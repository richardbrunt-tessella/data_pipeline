import unittest

from mrtarget.common.DataStructure import JSONSerializable
from mrtarget.common.ElasticsearchLoader import Loader

class SerializeStub(JSONSerializable):

    def __init__(self,
                 id,
                 **kwargs):
        self.id = id
        self.__dict__.update(**kwargs)

class ElasticsearchLoaderTestCase(unittest.TestCase):

    def test_put_string(self):
        loader = Loader(dry_run=True)
        loader.put('dummy-index',
                   'dummy-doctype',
                   'id',
                   '{"hello":"test"}')
        self.assertEquals(len(loader.cache),1)
        loader.flush()
        self.assertEquals(len(loader.cache),0)

    def test_put_jsonserializable(self):
        loader = Loader(dry_run=True)
        loader.put('dummy-index',
                   'dummy-doctype',
                   'id',
                   SerializeStub('test'))
        self.assertEquals(len(loader.cache),1)
        self.assertTrue(isinstance(loader.cache[0]['_source'], (str, unicode)))
        loader.flush()
        self.assertEquals(len(loader.cache),0)

    def test_many_put(self):
        loader = Loader(dry_run=True,
                        chunk_size=100)
        for i in range(150):
            loader.put('dummy-index',
                       'dummy-doctype',
                       'id',
                       '{"hello":"test"}')
        self.assertEquals(len(loader.cache),50)
        loader.close()
        self.assertEquals(len(loader.cache),0)



