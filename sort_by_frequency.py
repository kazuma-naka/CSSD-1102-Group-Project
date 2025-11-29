from collections import Counter


def sort_by_frequency(cids: list[int]) -> list[int]:
    """
    Sort the cids by frequency by descending order and remove duplicate items.

    NOTE: Your implementation may result in a list with different orders as in the example below as you are free to handle ties in any way you want
    In the autograder on PrairieLearn, we will NOT include examples with ties.
    """
    # If the input list is empty, simply return an empty list
    if not cids:
        return []

    # Count how many times each CID appears in the list
    counts = Counter(cids)

    # Record the first index where each CID appears.
    # This is used to break ties: if two CIDs have the same frequency,
    # the one that appears earlier in the original list comes first.
    first_index: dict[int, int] = {}
    for i, cid in enumerate(cids):
        if cid not in first_index:
            first_index[cid] = i

    # Extract unique CIDs (keys of the Counter)
    unique_cids = list(counts.keys())

    # Sort unique CIDs by:
    #   1) frequency in descending order  (-counts[cid])
    #   2) first occurrence index in ascending order  (first_index[cid])
    # This guarantees that:
    #   - higher frequency comes first
    #   - ties are resolved in a deterministic way
    unique_cids.sort(key=lambda cid: (-counts[cid], first_index[cid]))

    return unique_cids
