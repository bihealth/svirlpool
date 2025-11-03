import sqlite3

# Paths to your databases
db_A_path = "path/to/database_A.db"
db_B_path = "path/to/database_B.db"

# Connect to database B
conn_B = sqlite3.connect(db_B_path)
cur_B = conn_B.cursor()

# Attach database A to the connection of database B
cur_B.execute(f"ATTACH DATABASE '{db_A_path}' AS db_A")

# Copy the candidate_regions table from A to B
# Step 1: Check if candidate_regions table already exists in B
cur_B.execute(
    "SELECT name FROM sqlite_master WHERE type='table' AND name='candidate_regions'"
)
if not cur_B.fetchone():
    # Step 2: Create the table in B if it does not exist
    cur_B.execute(
        """
    CREATE TABLE candidate_regions AS SELECT * FROM db_A.candidate_regions
    """
    )
else:
    # If table exists, insert data
    cur_B.execute(
        """
    INSERT INTO candidate_regions (crID, candidate_region)
    SELECT crID, candidate_region FROM db_A.candidate_regions
    """
    )

# Detach database A from B
cur_B.execute("DETACH DATABASE db_A")

# Commit changes and close connection
conn_B.commit()
conn_B.close()

print("Table 'candidate_regions' has been copied from A to B successfully.")
