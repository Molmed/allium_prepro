import yaml
import os


class SubtypeThesaurus():
    def __init__(self):
        self._dict = {}
        self._allium_subtypes = []
        self._allium_groups = []
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
                if 'parent_id' in subtype:
                    self._allium_groups.append(subtype['parent_id'])
                self._dict[subtype['id']] = subtype['id']
                if 'aliases' in subtype:
                    for alias in subtype['aliases']:
                        self._dict[alias] = subtype['id']
            # Dedupe groups
            self._allium_groups = list(set(self._allium_groups))

    def is_allium_subtype(self, s):
        return s in self._allium_subtypes

    def allium_subtypes(self, include_groups=False):
        if include_groups:
            return self._allium_subtypes
        # Return subtypes that are not also in groups
        return list(set(self._allium_subtypes) - set(self._allium_groups))

    def translate(self, s):
        s = s.strip()
        return self._dict.get(s, f"UNRECOGNIZED_SUBTYPE: {s}")

    def translate_subtype_list(self, subtypes):
        translated = {}
        for entry in subtypes:
            split_subtypes = entry.split(',')
            if len(split_subtypes) > 1:
                translated_split_subtypes = []
                for x in split_subtypes:
                    translated_split_subtypes.append(self.translate(x))
                translated[entry] = ','.join(translated_split_subtypes)
            else:
                translated[entry] = self.translate(entry)
        return translated

    def translate_subtype_column(self, column):
        return column.replace(self.translate_subtype_list(column))

    def thesaurus(self):
        return self._dict
