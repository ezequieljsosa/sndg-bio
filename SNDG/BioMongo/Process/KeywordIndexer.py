"""
Created on Oct 2, 2015

@author: eze

http://norm.al/2009/04/14/list-of-english-stop-words/
"""

import re

STOP_WORDS = ["null", "none", "type", "a", "about", "above", "above", "across", "after", "afterwards", "again",
              "against", "all", "almost", "alone", "along", "already",
              "also", "although", "always", "am", "among", "amongst", "amoungst", "amount", "an", "and", "another",
              "any", "anyhow", "anyone", "anything", "anyway", "anywhere",
              "are", "around", "as", "at", "back", "be", "became", "because", "become", "becomes", "becoming", "been",
              "before", "beforehand", "behind", "being", "below", "beside",
              "besides", "between", "beyond", "bill", "both", "bottom", "but", "by", "call", "can", "cannot", "cant",
              "co", "con", "could", "couldnt", "cry", "de", "describe", "detail",
              "do", "done", "down", "due", "during", "each", "eg", "eight", "either", "eleven", "else", "elsewhere",
              "empty", "enough", "etc", "even", "ever", "every", "everyone",
              "everything", "everywhere", "except", "few", "fifteen", "fify", "fill", "find", "fire", "first", "five",
              "for", "former", "formerly", "forty", "found", "four", "from",
              "front", "full", "further", "get", "give", "go", "had", "has", "hasnt", "have", "he", "hence", "her",
              "here", "hereafter", "hereby", "herein", "hereupon", "hers",
              "herself", "him", "himself", "his", "how", "however", "hundred", "ie", "if", "in", "inc", "indeed",
              "interest", "into", "is", "it", "its", "itself", "keep", "last",
              "latter", "latterly", "least", "less", "ltd", "made", "many", "may", "me", "meanwhile", "might", "mill",
              "mine", "more", "moreover", "most", "mostly", "move", "much",
              "must", "my", "myself", "name", "namely", "neither", "never", "nevertheless", "next", "nine", "no",
              "nobody", "none", "noone", "nor", "not", "nothing", "now", "nowhere",
              "of", "off", "often", "on", "once", "one", "only", "onto", "or", "other", "others", "otherwise", "our",
              "ours", "ourselves", "out", "over", "own", "part", "per", "perhaps",
              "please", "put", "rather", "re", "same", "see", "seem", "seemed", "seeming", "seems", "serious",
              "several", "she", "should", "show", "side", "since", "sincere", "six",
              "sixty", "so", "some", "somehow", "someone", "something", "sometime", "sometimes", "somewhere", "still",
              "such", "system", "take", "ten", "than", "that", "the", "their",
              "them", "themselves", "then", "thence", "there", "thereafter", "thereby", "therefore", "therein",
              "thereupon", "these", "they", "thickv", "thin", "third", "this", "those",
              "though", "three", "through", "throughout", "thru", "thus", "to", "together", "too", "top", "toward",
              "towards", "twelve", "twenty", "two", "un", "under", "until", "up",
              "upon", "us", "very", "via", "was", "we", "well", "were", "what", "whatever", "when", "whence",
              "whenever", "where", "whereafter", "whereas", "whereby", "wherein",
              "whereupon", "wherever", "whether", "which", "while", "whither", "who", "whoever", "whole", "whom",
              "whose", "why", "will", "with", "within", "without", "would", "yet", "you",
              "your", "yours", "yourself", "yourselves", "the", "non"]

ALLOWED_SHORT_WORDS = ["Na", "Mg", "Al", "K", "Ca", "Cr", "Mn", "Fe", "Cu", "Zn", "Br"]
DEFAULT_IRRELEVANT_WORDS = ["according"]


class KeywordIndexer:

    def __init__(self):

        self.min_word_len = 3
        self.not_indexable_words = DEFAULT_IRRELEVANT_WORDS + STOP_WORDS
        self.allowed_short_words = ALLOWED_SHORT_WORDS

    def _extract_keywords(self, string):
        result = []
        for word in [x.lower() for x in re.findall('\w+', string.strip()) if
                     len(x) >= self.min_word_len or x in self.allowed_short_words]:
            if (word not in self.not_indexable_words) or (word in self.allowed_short_words):
                result.append(word)
        return result

    def extract_keywords(self, string):
        """
        Extracts a list of valid words from a string
        
        Parameters
        ----------
        string: containing creo or more words
        
        Returns
        -------
        list of "valid" words
        """
        if isinstance(string, list) or isinstance(string, tuple):
            result = []
            for stringx in string:
                result = result + self._extract_keywords(stringx)
        else:
            result = self._extract_keywords(string)
        return list(set(result))
