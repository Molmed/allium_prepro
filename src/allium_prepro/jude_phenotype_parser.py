import pandas as pd
import re
from .subtype_thesaurus import SubtypeThesaurus


class JudePhenotypeParser:
    # Definitions
    SUBTYPE_PRIMARY = 1
    SUBTYPE_SECONDARY = 2
    SUBTYPE_LEVEL_DELIMITER = ','

    def __init__(self,
                 prefix,
                 pheno_tsv,
                 output_dir):
        self._prefix = prefix
        self._pheno_tsv = pheno_tsv
        self._output_dir = output_dir
        self._output_file_path = f'{output_dir}/{prefix}.pheno.allium.csv'
        self._subtype_thesaurus = SubtypeThesaurus()

        # Track unknown subtypes
        self.unknown_primary_subtypes = {}
        self.unknown_secondary_subtypes = {}

        # Keep the underlying dataframe
        self.df = None

    @staticmethod
    def split_attr_diagnosis(s):
        # Define the regular expression pattern
        pattern = r'Lineage:(.*?),Primary_subtype:(.*?)(?:,Secondary_subtype:(.*))?$'

        # Use re.match() to apply the pattern
        match = re.match(pattern, s)

        # Extract the groups from the match object
        lineage = match.group(1)
        primary_subtype = match.group(2)
        secondary_subtype = match.group(3) \
            if match.group(3) is not None else ''

        return lineage, primary_subtype, secondary_subtype

    def alliumify_subtype(self, s, level=SUBTYPE_PRIMARY):
        # If level is SUBTYPE_PRIMARY and s is empty string, return "B-other"
        if level == JudePhenotypeParser.SUBTYPE_PRIMARY and s == '':
            return 'B-other'

        # Strip trailing ?s
        s = s.rstrip('?')

        # If subtype is recognized, return it
        if self._subtype_thesaurus.is_allium_subtype(s):
            return s

        # For secondary subtypes, return empty string if s is empty
        if s == '':
            return ''

        # Get translation
        translation = self._subtype_thesaurus.translate(s)

        if translation.startswith('UNRECOGNIZED_SUBTYPE'):
            if level == JudePhenotypeParser.SUBTYPE_PRIMARY:
                self.unknown_primary_subtypes[s] = \
                    self.unknown_primary_subtypes.get(s, 0) + 1
            elif level == JudePhenotypeParser.SUBTYPE_SECONDARY and s != '':
                self.unknown_secondary_subtypes[s] = \
                    self.unknown_secondary_subtypes.get(s, 0) + 1
            return ""  # Return empty string for unrecognized subtypes

        return translation

    def parse(self):
        print("Parsing St. Jude phenotype data...")

        df = pd.read_csv(self._pheno_tsv, delimiter='\t')

        # Drop AML cases
        df = df[df['attr_diagnosis'] != 'AML']

        # Drop unnecessary columns
        df = df[['sample_name', 'attr_diagnosis']]

        # Rename sample_name to id
        df = df.rename(columns={'sample_name': 'id'})

        for index, row in df.iterrows():
            lineage, primary_subtype, secondary_subtype = \
                JudePhenotypeParser.split_attr_diagnosis(row['attr_diagnosis'])

            # ALLIUM only has one T-ALL subtype
            if lineage == 'T':
                primary_subtype = 'T-ALL'

            # Standardize to ALLIUM conventions
            primary_subtype = self.alliumify_subtype(
                primary_subtype, level=self.SUBTYPE_PRIMARY)
            secondary_subtype = self.alliumify_subtype(
                secondary_subtype, level=self.SUBTYPE_SECONDARY)

            # Remove empty strings, sort and join
            final_subtype = self.SUBTYPE_LEVEL_DELIMITER.join(sorted([
                x for x in [primary_subtype, secondary_subtype] if x != '']))

            # Write final_subtype to 'subtype' column
            df.at[index, 'subtype'] = final_subtype

        # Drop 'attr_diagnosis' column
        df = df.drop(columns=['attr_diagnosis'])

        # Store the df
        self.df = df

        # Save summary
        self.save_summary()
        self.print_summary()

        self.df.to_csv(self._output_file_path, sep=';', index=False)

    def save_summary(self):
        with open(f'{self._output_dir}/{self._prefix}.pheno_summary.txt', 'w') as f:
            f.write("Unknown primary subtypes: %s\n" %
                    self.unknown_primary_subtypes)
            f.write("Total cases with unknown primary subtypes: %d\n" %
                    sum(self.unknown_primary_subtypes.values()))
            f.write("Unknown secondary subtypes: %s\n" %
                    self.unknown_secondary_subtypes)
            f.write("Total cases with unknown secondary subtypes: %d\n" %
                    sum(self.unknown_secondary_subtypes.values()))
            f.write("Number of cases remaining: %d\n" % self.df.shape[0])

    def print_summary(self):
        with open(f'{self._output_dir}/{self._prefix}.pheno_summary.txt', 'r') as f:
            print(f.read())
