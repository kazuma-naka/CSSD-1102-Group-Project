import csv
from io import StringIO


def remove_non_drug(sorted_cids: list[int], save_result: bool = False, mol_name: str = "") -> list[int]:
    """
    Given a list of sorted cids, remove the non-drug structures based on Lipinski's rule of five
    The queries involved in this function will be time-consuming, you will need to create mulitple chunks.
    Ignore the entries that don't have all five properties available.

    Args:
         sorted_cids - a list of sorted cids
         save_result - whether to save the result to a csv file. The default is False. Please make sure this is set to the False when send your code to the autograder
         mol_name    - the name of the compound. This will be used as the name of the csv file generated. You must provide it when save_result is True
    """
    if not sorted_cids:
        return []

    # Lipinski-like conditions
    MAX_MW = 500.0
    MAX_XLOGP = 5.0
    MAX_HBD = 5
    MAX_HBA = 10

    # This dict will store properties for each CID that passes Lipinski rules
    # so that we can optionally write them to a CSV file.
    druglike_props: dict[int, tuple[float, float, int, int]] = {}

    # Split the CID list into chunks to avoid overly large single requests.
    chunks = get_chunks(sorted_cids, chunk_size=100)

    for chunk in chunks:
        if not chunk:
            continue

        # Build a comma-separated list of CIDs for this chunk
        cid_str = ",".join(str(cid) for cid in chunk)

        # Prepare the property list to query from PubChem
        properties = "MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount"

        # Example URL:
        # https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/962,910/property/MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount/csv
        url = generate_url(
            f"compound/cid/{cid_str}",
            f"property/{properties}",
            "csv",
        )

        # Query PubChem and get the CSV text
        csv_text = make_query(url)

        # Parse the CSV text using csv.DictReader for convenience
        reader = csv.DictReader(StringIO(csv_text))

        for row in reader:
            # Skip rows that are missing any of the required columns
            try:
                cid_value = row.get("CID", "").strip()
                mw_value = row.get("MolecularWeight", "").strip()
                xlogp_value = row.get("XLogP", "").strip()
                hbd_value = row.get("HBondDonorCount", "").strip()
                hba_value = row.get("HBondAcceptorCount", "").strip()

                # If any field is empty or missing, ignore this entry
                if not (cid_value and mw_value and xlogp_value and hbd_value and hba_value):
                    continue

                cid = int(cid_value)
                mw = float(mw_value)
                xlogp = float(xlogp_value)
                # Sometimes these may be represented as floats in the CSV, so cast via float -> int
                hbd = int(float(hbd_value))
                hba = int(float(hba_value))
            except (ValueError, KeyError):
                # If parsing fails for any reason, skip this row
                continue

            # Apply Lipinski's rule of five
            if (
                mw <= MAX_MW
                and xlogp <= MAX_XLOGP
                and hbd <= MAX_HBD
                and hba <= MAX_HBA
            ):
                druglike_props[cid] = (mw, xlogp, hbd, hba)

    # Preserve the original order from sorted_cids,
    # but only keep CIDs that passed the Lipinski filter.
    druglike_cids = [cid for cid in sorted_cids if cid in druglike_props]

    # Optionally write results to a CSV file
    if save_result:
        if not mol_name:
            # A name is required when saving to file, to determine the file name
            raise ValueError(
                "mol_name must be provided when save_result is True")

        filename = f"{mol_name}.csv"
        with open(filename, "w", newline="") as f:
            writer = csv.writer(f)
            # Header row
            writer.writerow(
                ["CID", "MolecularWeight", "XLogP",
                    "HBondDonorCount", "HBondAcceptorCount"]
            )
            # Write properties in the same order as druglike_cids
            for cid in druglike_cids:
                mw, xlogp, hbd, hba = druglike_props[cid]
                writer.writerow([cid, mw, xlogp, hbd, hba])

    return druglike_cids
