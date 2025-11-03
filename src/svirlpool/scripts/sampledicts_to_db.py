import argparse
import json
import sqlite3
from pathlib import Path

from . import util


def add_sampledicts_to_db(database: Path, sampledicts: Path):
    sd = util.load_sampledicts(sampledicts)
    sd_txt = json.dumps(sd)
    # write sd_txt to database
    with sqlite3.connect(str(database)) as conn:
        conn.execute("DROP TABLE IF EXISTS sampledicts")
        conn.execute(
            """CREATE TABLE sampledicts (
                        sampledicts TEXT)"""
        )
        c = conn.cursor()
        c.execute("INSERT INTO sampledicts (sampledicts) VALUES (?)", (sd_txt,))
        c.close()
        conn.commit()


def read_sampledicts_from_db(database: Path) -> dict:
    try:
        with sqlite3.connect(str(database)) as conn:
            c = conn.cursor()
            c.execute("SELECT sampledicts FROM sampledicts")
            result = c.fetchone()
            conn.commit()
            c.close()
    except:
        raise ValueError(f"could not read sampledicts from database {database}")
    return json.loads(result[0])


def run(args, **kwargs):
    add_sampledicts_to_db(database=args.output, sampledicts=args.sampledicts)


def get_parser():
    parser = argparse.ArgumentParser(
        description="Call SVs from the final database and write them to a vcf file."
    )
    parser.add_argument(
        "-i",
        "--sampledicts",
        type=Path,
        required=True,
        help="Path to the sampledicts file.",
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True, help="Path to the output vcf file."
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
