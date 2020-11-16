from mongoengine import EmbeddedDocument
from mongoengine.fields import StringField, ListField, DynamicField


class SeqColDruggabilityParamTypes():
    number = "number"
    value = "value"

    values = [number, value]


class   SeqColDruggabilityParam(EmbeddedDocument):
    '''
    '''

    overexpressed = ["stress", "starvation", "infection", "hypoxia"]
    default_params = [("essentiality",
                       "Critical for the organism survival (https://www.ncbi.nlm.nih.gov/pubmed/26791267)", "protein",
                       SeqColDruggabilityParamTypes.value, ["true", "false"],
                       "avg", "equal", "true"),
                      ("human_offtarget", """
                     This score reflects the results of a blastp search of the pathogen protein in the human proteome database (ncbi accession GCF_000001405.36)
                     with the scale 1 - max(alignment identity), so when a protein has no hit in the human proteome,  
                     the value is 1, and if it has 2 hits, one with an identity of 0.4 and other with 0.6, the score is 0.4 (human_offtarget = 1 - 0.6, uses the max identity).                     
                     """.strip(), "protein", SeqColDruggabilityParamTypes.number, None,
                       "max", "<", 0.4),
                      ("hit_in_deg", "Has a hit in Database of Essential Genes", "protein",
                       SeqColDruggabilityParamTypes.value, ["Yes", "No"],
                       "avg", "equal", "Yes")
                      ]
    '''
    name,description,target,_type,options,defaultGroupOperation,defaultOperation,defaultValue
    '''
    for cond in overexpressed:
        default_params.append(("overexpression_" + cond,
                               "Overexpressed in model of " + cond + " (https://www.ncbi.nlm.nih.gov/pubmed/26791267)",
                               "protein", SeqColDruggabilityParamTypes.value, ["true", "false"],
                               "avg", "equal", "true"))

    MAX_OPTIONS = 20

    meta = {'allow_inheritance': True, 'strict': False}
    name = StringField(required=True)
    description = StringField(default="")
    type = StringField(default="value", choices=["number","value"])
    target = StringField(default="protein")
    uploader = StringField(default="")
    options = ListField(DynamicField())
    _class = StringField(default="ar.com.bia.entity.SeqCollectionDoc")
    defaultGroupOperation = StringField(required=False)
    defaultOperation = StringField(required=False)
    defaultValue = DynamicField(required=False)

    def isValid(self):
        return self.type == SeqColDruggabilityParamTypes.number.value or len(
            self.options) <= SeqColDruggabilityParam.MAX_OPTIONS

    def __str__(self):
        return "%s type='%s' target='%s'" % (self.name, self.type, self.target)
