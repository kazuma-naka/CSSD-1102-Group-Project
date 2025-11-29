import csv
from io import StringIO


def sort_and_retrieve_property(
    cids: list[int],
    compound_property: str,
    properties_to_include: list[str],
    descending: bool = True,
):
    """
    Given a list of cids, retrieve a few properties (e.g., properties_to_include)
    and sort them based on a particular compound_property.

    If descending is True, sort them in descending order. Otherwise, sort them in ascending order.
    The queries involved in this function will be time-consuming; multiple chunks are used.

    This function is mainly used to support CLO5 questions (e.g., sorting by MolecularWeight,
    then checking XLogP among the top 10 entries).

    Returns:
        A list of dictionaries, where each dictionary corresponds to one row from PubChem
        (keys are column names such as "CID", "MolecularWeight", "XLogP", etc.).
    """
    # If there is nothing to query or no properties requested, return an empty list
    if not cids or not properties_to_include:
        return []

    # Build the property list string for PubChem, e.g. "MolecularWeight,XLogP"
    properties_str = ",".join(properties_to_include)

    # This list will accumulate all rows (as dictionaries) returned from PubChem
    all_rows: list[dict[str, str]] = []

    # Use chunking to avoid sending too many CIDs in a single request
    chunks = get_chunks(cids, chunk_size=100)

    for chunk in chunks:
        if not chunk:
            continue

        # Build comma-separated CID list for this chunk
        cid_str = ",".join(str(cid) for cid in chunk)

        # Example URL:
        # https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/962,910/property/MolecularWeight,XLogP/csv
        url = generate_url(
            f"compound/cid/{cid_str}",
            f"property/{properties_str}",
            "csv",
        )

        # Query PubChem and get raw CSV text
        csv_text = make_query(url)

        # Use DictReader to parse CSV into dictionaries keyed by column name
        reader = csv.DictReader(StringIO(csv_text))

        for row in reader:
            # Ensure there is at least a CID field and the sorting property
            cid_value = row.get("CID", "").strip()
            sort_value = row.get(compound_property, "").strip()

            # Skip rows missing CID or the sorting property
            if not cid_value or not sort_value:
                continue

            # Keep the row as-is (string values). Conversion is done in the sort key.
            all_rows.append(row)

    # Define a helper function to extract a sortable key from each row
    def _sort_key(row: dict[str, str]):
        value_str = row.get(compound_property, "").strip()
        # Try to interpret the value as a float (typical for numeric properties)
        try:
            return float(value_str)
        except ValueError:
            # If conversion fails (e.g., non-numeric property), fall back to string
            return value_str

    # Sort the accumulated rows by the chosen compound_property
    all_rows.sort(key=_sort_key, reverse=descending)

    # Return the sorted data so the caller can:
    # - inspect it,
    # - print it,
    # - or save it to a CSV if needed.
    return all_rows
