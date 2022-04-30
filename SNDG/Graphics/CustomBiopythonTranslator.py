from dna_features_viewer import BiopythonTranslator


class FeatureTypeFilter:
    def __init__(self, feature_types):
        self.feature_types = feature_types

    def filter_feature(self, feature):
        if feature.type in self.feature_types:
            return True
        return False


class FeatureTypeColourer:
    def __init__(self, color_map, default_color=None):
        self.color_map = color_map
        self.default_color = default_color

    def feature_color(self, feature):
        return self.color_map[feature.type] if feature.type in self.color_map else (
            self.default_color if self.default_color else None)


class FeatureQualifierColourer:
    def __init__(self, qualifier, color_map, default_color=None):
        self.qualifier = qualifier
        self.color_map = color_map
        self.default_color = default_color

    def feature_color(self, feature):
        if self.qualifier in feature.qualifiers:
            qual_value = feature.qualifiers[self.qualifier][0]

            if qual_value in self.color_map:
                return self.color_map[qual_value]
            elif self.default_color:
                return self.default_color
        return None


class FeatureTypeQualifierLabeler:
    def __init__(self, type_qualifiers_map, default_label=None):
        self.type_qualifiers_map = type_qualifiers_map
        self.default_label = default_label

    def feature_label(self, feature):
        if feature.type in self.type_qualifiers_map:
            qualifiers = self.type_qualifiers_map[feature.type]
            for qualifier in qualifiers:
                if qualifier in feature.qualifiers and feature.qualifiers[qualifier][0].strip():
                    qualifier_value = feature.qualifiers[qualifier][0].strip()
                    return qualifier_value
        if self.default_label:
            return self.default_label

class FeatureQualifierValueLabeler:
    def __init__(self, qualifier, qualifier_label_map):
        self.qualifier = qualifier
        self.qualifier_label_map = qualifier_label_map

    def feature_label(self, feature):
        if self.qualifier in feature.qualifiers and feature.qualifiers[self.qualifier][0].strip():
            qualifier_value = feature.qualifiers[self.qualifier][0].strip()
            if qualifier_value in self.qualifier_label_map:
                return self.qualifier_label_map[qualifier_value]


class CustomBiopythonTranslator(BiopythonTranslator):


    def __init__(self, *args, **kwargs):
        super(BiopythonTranslator, self).__init__(*args, **kwargs)
        self.feature_filters = []
        self.feature_colourers = []
        self.feature_labelers = []




    def add_filter(self, feature_filter):
        self.feature_filters.append(feature_filter)

    def add_colourer(self, colourer):
        self.feature_colourers.append(colourer)

    def add_labeler(self, labeler):
        self.feature_labelers.append(labeler)

    @staticmethod
    def buildDefault():
        cbt = CustomBiopythonTranslator()
        cbt.add_filter(FeatureTypeFilter(["CDS", "GAP", "REGION"]))
        cbt.add_colourer(FeatureTypeColourer({"CDS": "red", "GAP": "black", "REGION": "purple"}, default_color="grey"))
        cbt.add_labeler(FeatureTypeQualifierLabeler({"CDS": ["gene", "locus_tag"], "REGION": ["name"]}))
        return cbt

    def compute_feature_color(self, feature):
        for colourer in self.feature_colourers:
            color = colourer.feature_color(feature)
            if color:
                return color
        return "white"

    def compute_feature_label(self, feature):
        for labeler in self.feature_labelers:
            label = labeler.feature_label(feature)
            if label:
                return label
        return None

    def compute_filtered_features(self, features):
        for f in features:
            for feature_filter in self.feature_filters:
                if feature_filter.filter_feature(f):
                    yield f
