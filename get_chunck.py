def get_chunks(cids: list[int], chunk_size: int = 100) -> list[list[int]]:
    """
    Divide the list into multiple chunks of certain chunk_size.
    """
    # chunk_size must be a positive integer
    if chunk_size <= 0:
        raise ValueError("chunk_size must be a positive integer")

    chunks: list[list[int]] = []

    # Iterate over the list in steps of chunk_size
    # For example, if chunk_size = 3 and len(cids) = 8:
    #   i = 0 -> cids[0:3]
    #   i = 3 -> cids[3:6]
    #   i = 6 -> cids[6:9] (last slice will be shorter if needed)
    for i in range(0, len(cids), chunk_size):
        # Append a slice of length at most chunk_size
        chunks.append(cids[i:i + chunk_size])

    return chunks
