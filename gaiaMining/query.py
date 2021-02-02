# Functions to query Gaia EDR3
# =========================================================================
# Imports
# =========================================================================
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia


if __name__ == "__main__":
  # Perform some tests to ensure basic functionalities are present
  # Options for printing, etc
  Gaia.ROW_LIMIT = 8

  # Basic query through coordinates, width and height
  coord = SkyCoord(ra=280, dec=60, unit=(u.degree, u.degree), frame='icrs')
  width = u.Quantity(0.1, u.deg)
  height = u.Quantity(0.1, u.deg)
  print("=================================================================")
  print(f"Basic query with coord {coord}, width {width} and height {height}")
  print("=================================================================")
  r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
  r.pprint()

  # Cone search centred at specified RA/Dec coords with provided radius
  # Using 'coord' as defined above
  radius = u.Quantity(1.0, u.deg)
  j = Gaia.cone_search_async(coord, radius)
  print("=================================================================")
  print(f"Cone search with coord {coord} and radius {radius}")
  print("=================================================================")
  r = j.get_results()
  r.pprint()

  # Getting public tables metadata
  print("=================================================================")
  print("Table names metadata")
  print("=================================================================")
  tables = Gaia.load_tables(only_names=True)
  filter_string = 'gaiaedr3'
  filtered_tables = []

  for table in tables:
    # Filter by 'gaiadr3'
    name = table.get_qualified_name()
    if filter_string in name:
      print(name)
      filtered_tables.append(table)

  # Loading single data and accessing column names
  print("=================================================================")
  print("Single table columns")
  print("=================================================================")
  # Just a demo of loading a table from its name
  table_query_name = '.'.join(filtered_tables[1].get_qualified_name().split('.')[1:])
  gaiaedr3_table = Gaia.load_table(table_query_name)

  # Query and save results to file
  # !!!
  # Dumping to file results in error
  # AttributeError: 'str' object has no attribute 'read'
  
  test_out = './out/gaiaedr3_test.csv'
  print("=================================================================")
  print(f"Single query and saving results to {test_out}")
  print("=================================================================")
  query = "select top 100 "\
    "solution_id,ref_epoch,ra_dec_corr,astrometric_n_obs_al, "\
    "matched_observations,duplicated_source,phot_variable_flag "\
    f"from {table_query_name} by source_id"
  job = Gaia.launch_job(query, dump_to_file=True)
  # !!!
  