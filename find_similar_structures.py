import doctest
from generate_url import generate_url, make_query


def find_similar_structures(cid: int) -> list[int]:
    """
    Given a valid cid, return a list of similar structures based on fastsimilarity_2d

    >>> find_similar_structures(962)
    [962, 961, 24602, 104752, 22247451, 10129877, 157350, 6914120, 123332, 105142, 137349153, 139859, 5460554, 10197601, 16048639, 6914119, 17985087, 21871715, 22146520, 44144404, 57470146, 129631210, 171528, 10214376, 20266850, 57422081, 72499223, 5460604, 12676391, 21338204, 23560202, 57418915, 57634601, 57937911, 58422532, 58609138, 58750500, 58828612, 59045931, 59049839, 59102850, 59367517, 59984341, 78060557, 85595863, 85602483, 139617572, 157621312, 158956564, 160528969, 161115002, 163575618, 139505, 12161503, 13644079, 15593902, 17769776, 19818716, 20037327, 71327550, 71333220, 72700706, 76071083]
    """
    # Build URL for fast 2D similarity search
    url = generate_url(
        f"compound/fastsimilarity_2d/cid/{cid}",
        "cids",
        "txt",
    )
    # Query PubChem
    text = make_query(url)
    # Split into lines, drop empty strings, convert to int
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    return [int(line) for line in lines]


doctest.testmod()
