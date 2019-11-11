import sqlite3, time
import pandas as pd

'''
Results with ~40000 entries, search for "Coils":

Length of SQLhits: 4557
SQL performance: 0.053612709045410156
Length of Pandashits: 4557
Pandas performance: 15.488629341125488
'''

def benchmark_functional(df, filename="csvdescriptions.sql"):
	df = pd.DataFrame(df).T
	conn = sqlite3.connect(filename)

	df.to_sql("functional", conn, if_exists='replace', index=True)

	# Measurements
	query = "Coil"

	# SQL
	start = time.time()
	# String generator-function needs to be added to cope with dynamic presence of df-headers (TIGRFAM LIKE '%{0}%'.... etc)
	sqldf = pd.read_sql_query("select * from functional where Coils LIKE '%{0}%' or SUPERFAMILY LIKE '%{0}%' or TIGRFAM LIKE '%{0}%' or Hamap LIKE '%{0}%' or SMART LIKE '%{0}%' or ProSitePatterns LIKE '%{0}%' or PIRSF LIKE '%{0}%';".format(query), conn)
	print("Length of SQLhits: {}".format(len(sqldf)))
	end = time.time()
	print("SQL performance: {}".format(end - start))

	# Pandas
	start = time.time()
	pandasdf = df.where(df.apply(lambda row: row.astype(str).str.contains(query, case=False).any(), axis=1, result_type='broadcast')).dropna(how='all')
	print("Length of Pandashits: {}".format(len(pandasdf)))
	end = time.time()
	print("Pandas performance: {}".format(end - start))
