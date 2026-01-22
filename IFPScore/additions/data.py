"""
Additional functionality.
"""
from typing import Literal

import pandas as pd
from qsprpred.data.data import MoleculeTable, DescriptorTable
from rdkit import Chem
from qsprpred.data.sources.papyrus import Papyrus
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from spock.storage.base import ChemStore
from tqdm.auto import tqdm

# different grid visualizations
standard_grid = Chem.Draw.MolsToGridImage
def interactive_grid(mols, *args, molsPerRow=5, **kwargs):
    """
    install mols2grid with pip to use
    """

    import mols2grid

    return mols2grid.display(mols, *args, n_cols=molsPerRow, **kwargs)


def smiles_to_grid(smiles, *args, mols_per_row=5, **kwargs):
    mols = []
    for smile in smiles:
        try:
            m = Chem.MolFromSmiles(smile)
            if m:
                AllChem.Compute2DCoords(m)
                mols.append(m)
            else:
                raise Exception(f'Molecule empty for SMILES: {smile}')
        except Exception as exp:
            pass

    return interactive_grid(mols, *args, molsPerRow=mols_per_row, **kwargs)


class MolTable(MoleculeTable):

    @staticmethod
    def matchMolToSMARTS(
            mol: Chem.Mol | str,
            smarts: list[str],
            operator: Literal["or", "and"] = "or",
            use_chirality: bool = False
    ):
        """
        Check if a molecule matches a SMARTS pattern.

        Args:
            mol: Molecule to check.
            smarts: SMARTS pattern to check.
            operator: Whether to use an "or" or "and" operator on patterns. Defaults to "or".
            use_chirality: Whether to use chirality in the search.

        Returns:
            (bool): True if the molecule matches the pattern, False otherwise.
        """
        mol = Chem.MolFromSmiles(mol) if isinstance(mol, str) else mol
        ret = False
        for smart in smarts:
            ret = mol.HasSubstructMatch(Chem.MolFromSmarts(smart), useChirality=use_chirality)
            if operator == "or":
                if ret:
                    return True
            elif operator == "and":
                if ret:
                    ret = True
                else:
                    return False
        return ret

    def toGrid(self, mols_per_row=5):
        return smiles_to_grid(self.df[self.smilesCol], mols_per_row=mols_per_row)

    def searchWithIndex(self, index: pd.Index, name: str | None = None):
        """
        Create a new table from a list of indices.

        Args:
            index(pd.Index): Indices in this table to create the new table from.
            name(str): Name of the new table. Defaults to the name of the old table, plus the `_new` suffix.

        Returns:
            MolTable: A new table with the molecules from the old table with the given indices.
        """
        name = f"{self.name}_new" if name is None else name
        ret = MolTable(
            name=name,
            df=self.df.loc[index, :],
            smiles_col=self.smilesCol,
            add_rdkit=False,
            store_dir=self.storeDir,
            overwrite=True,
            n_jobs=self.nJobs,
            chunk_size=self.chunkSize,
            drop_invalids=False,
            index_cols=self.indexCols
        )
        for table, calc in zip(self.descriptors, self.descriptorCalculators):
            ret.descriptors.append(
                DescriptorTable(
                    calc,
                    name_prefix=name,
                    df=table.getDF().loc[index, :],
                    store_dir=table.storeDir,
                    overwrite=True,
                    key_cols=table.indexCols,
                    n_jobs=table.nJobs,
                    chunk_size=table.chunkSize
                )
            )
        ret.descriptorCalculators = self.descriptorCalculators
        return ret

    def searchOnColumn(self, col_name: str, values: list[str], name: str | None = None, exact=False):
        mask = [False] * len(self.df)
        for value in values:
            mask = mask | (self.df[col_name].str.contains(value)) if not exact else mask | (self.df[col_name] == value)
        matches = self.df.index[mask]
        return self.searchWithIndex(matches, name)

    def searchWithSMARTS(
            self,
            patterns: list[str],
            operator: Literal["or", "and"] = "or",
            use_chirality: bool = False,
            name: str | None = None
    ):
        """
        Search the molecules in the table with a SMARTS pattern.

        Args:
            patterns: List of SMARTS patterns to search with.
            operator (object): Whether to use an "or" or "and" operator on patterns. Defaults to "or".
            use_chirality: Whether to use chirality in the search.
            name: Name of the new table. Defaults to the name of the old table, plus the `smarts_searched` suffix.

        Returns:
            (MolTable): A dataframe with the molecules that match the pattern.
        """
        matches = self.df.index[self.df[self.smilesCol].apply(
            lambda x: self.matchMolToSMARTS(x, patterns, operator=operator, use_chirality=use_chirality)
        )]

        return self.searchWithIndex(matches, name=f"{self.name}_smarts_searched" if name is None else name)

    def getSummary(self):
        """
        Make a summary of some statistics from the table

        Returns:
            (pd.DataFrame): A dataframe with the summary statistics.

        """
        # create a summary dictionary
        summary = {}
        # count number of compounds per target
        summary["mols_per_target"] = self.df.groupby("accession").count()["InChIKey"].to_dict()
        # count number of unique molecules
        summary["mols_per_target_unique"] = self.df.groupby("accession").aggregate(lambda x: len(set(x)))[
            "InChIKey"].to_dict()
        return pd.DataFrame(summary)


class PapyrusExtended(Papyrus):

    def getData(self, *args, name=None, output_dir=None, **kwargs):
        try:
            return MolTable(name=name, store_dir=output_dir or self.dataDir)
        except Exception as exp:
            pass
        ret = super().getData(*args, name=name, output_dir=output_dir, **kwargs)
        ret.__class__ = MolTable
        return ret


def moltable_to_store(table: MoleculeTable, storage: ChemStore, add_cols=None, enum_stereo=True):
    for idx, row in tqdm(table.df.iterrows(), total=len(table.df), desc="Adding molecules to Spock storage"):
        metadata = row[add_cols].to_dict()
        for col in table.indexCols:
            metadata[col] = row[col]
        metadata["Index"] = table.indexCols
        if enum_stereo:
            stereoisomers = tuple(EnumerateStereoisomers(Chem.MolFromSmiles(row[table.smilesCol])))
            for stereoisomer in stereoisomers:
                storage.add_mol_from_smiles(Chem.MolToSmiles(stereoisomer, isomericSmiles=True), metadata=metadata, update_existing="append")
        else:
            storage.add_mol_from_smiles(row[table.smilesCol], metadata=metadata, update_existing="append")
