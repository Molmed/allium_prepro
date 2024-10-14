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

    def translate(self, s):
        return self._dict.get(s, None)
