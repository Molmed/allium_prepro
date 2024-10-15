import yaml
import os


class SubtypeThesaurus():
    def __init__(self):
        self._dict = {}
        self._allium_subtypes = []
        # Get current path
        current_path = os.path.dirname(os.path.realpath(__file__))
        CONSTANTS_YML = f'{current_path}/subtypes.yml'
        # Load subtype mapping
        with open(CONSTANTS_YML, 'r') as f:
            constants = yaml.safe_load(f)
            subtypes = constants['subtypes']
            # Loop through subtypes
            for subtype in subtypes:
                self._allium_subtypes.append(subtype['id'])
                if 'aliases' in subtype:
                    for alias in subtype['aliases']:
                        self._dict[alias] = subtype['id']

    def is_allium_subtype(self, s):
        return s in self._allium_subtypes

    def translate(self, s, return_self_as_default=False):
        s = s.strip()
        if return_self_as_default:
            return self._dict.get(s, s)
        translated = self._dict.get(s, f"UNRECOGNIZED_SUBTYPE: {s}")
        return translated

    def translate_subtype_list(self, subtypes):
        translated = {}
        for entry in subtypes:
            split_subtypes = entry.split(',')
            if len(split_subtypes) > 1:
                translated_split_subtypes = []
                for x in split_subtypes:
                    translated_split_subtypes.append(self.translate(x, True))
                translated[entry] = ','.join(translated_split_subtypes)
            else:
                translated[entry] = self.translate(entry, True)
        return translated

    def thesaurus(self):
        return self._dict
