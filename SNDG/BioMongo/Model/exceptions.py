'''
Created on Jun 15, 2016

@author: eze
'''

class NotFoundException(Exception):
    '''
    classdocs
    '''


    def __init__(self, element):
        '''
        Constructor
        '''
        self.elementNotFound = element
        
    def __str__(self, *args, **kwargs):
        return "NotFoundException(%s)" % self.elementNotFound
        
        