'''
Created on Jul 1, 2011

@author: tel
'''         
class Data(object):
    def __init__(self, data_tuples=None, dtype = None, **kwargs):
        if data_tuples == None:
            data_tuples = []
        if dtype == None:
            dtype = 'Base'
        self.data_list = []
        self.dtype = dtype
        self.data = []
        self.AddDataTuple(data_tuples)
        for key in kwargs:
            self.__setattr__(key, kwargs[key])
    
    def Get(self, key):
        return object.__getattribute__(self, key)
    
    def AddDataTuple(self, data_tuples):
        for tup in data_tuples:
            self.data_list.append(tup[0])
            self.__setattr__(tup[0], tup[1])
    
    def GetDataTuples(self, exclude=(), include=()):
        temp = []
        if exclude==() and include==():
            for key in self.data_list:
                temp.append((key, object.__getattribute__(self, key)))
        elif exclude!=() and include==():
            for key in self.data_list:
                if key not in exclude:
                    temp.append((key, object.__getattribute__(self, key)))
        elif exclude==() and include!=():
            for key in self.data_list:
                if key in include:
                    temp.append((key, object.__getattribute__(self, key)))
        else:
            for key in self.data_list:
                if key in include:
                    if key not in exclude:
                        temp.append((key, object.__getattribute__(self, key)))
        return temp

    def GenDataTuples(self, exclude=(), include=()):
        if exclude==() and include==():
            for key in self.data_list:
                yield (key, object.__getattribute__(self, key))
        elif exclude!=() and include==():
            for key in self.data_list:
                if key not in exclude:
                    yield (key, object.__getattribute__(self, key))
        elif exclude==() and include!=():
            for key in self.data_list:
                if key in include:
                    yield (key, object.__getattribute__(self, key))
        else:
            for key in self.data_list:
                if key in include:
                    if key not in exclude:
                        yield (key, object.__getattribute__(self, key))
                
    def RemoveDataTuple(self, name):
        self.__delattr__(name)
        self.data_list.pop(self.data_list.index(name))

    def PopDataTuple(self, name):
        val = object.__getattribute__(self, name)
        self.RemoveDataTuple(name)
        return val

    #the data attribute and AddData method allow for data to contain data in a hierarchy
    #important for child class behavior
    def AddData(self, data):
        self.data.append(data)
