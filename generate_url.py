import doctest
import requests


def make_query(url: str) -> str:
    """
    Given a URL, make a requst to the ChemPub online database
    The interpreter complains if the escape character is not used (this is add only for the purpose of accommodating doctest)

    >>> make_query("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/water/cids/txt")
    '962\\n'
    """
    result = requests.get(url)

    # You can find what the code represents here: https://pubchem.ncbi.nlm.nih.gov//docs/pug-rest#section=Status-Codes
    assert result.status_code == 200, f"Error with url {url}\nstatus_code: {result.status_code}\nTo interpret the status code, please visit https://pubchem.ncbi.nlm.nih.gov//docs/pug-rest#section=Status-Codes\ncontent: {result.content}"

    return result.text


def generate_url(input: str, operation: str, output: str, options: str = "") -> str:
    """
    Generates a URL for the PubChem PUG REST API based on the provided input, operation, output format, and options.

    >>> generate_url("compound/name/butadiene", "property/HeavyAtomCount", "txt")
    'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/butadiene/property/HeavyAtomCount/txt'

    >>> generate_url("compound/cid/129825914,129742624,129783988", "sids", "json", "list_return=grouped")
    'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/129825914,129742624,129783988/sids/json?list_return=grouped'
    """
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
    url = f"{base}{input}/{operation}/{output}"
    if options:
        url = f"{url}?{options}"
    return url


def name_to_cid(name: str) -> int:
    """
    Given the name of a compound, return the cid

    >>> name_to_cid("water")
    962
    """
    # input: compound/name/{name}
    # operation: cids
    # output: txt
    url = generate_url(f"compound/name/{name}", "cids", "txt")
    text = make_query(url)
    cid_str = text.strip()
    return int(cid_str)


doctest.testmod()
