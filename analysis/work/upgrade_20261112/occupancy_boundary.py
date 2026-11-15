import sys
from pathlib import Path
sys.path.insert(0,'F:/analysis/climate_biliary_20261102/work/pylibs')
import duckdb
base=Path(__file__).resolve().parents[2]
out=base/'outputs/upgrade_20261112'
work=base/'work/upgrade_20261112'
db=duckdb.connect()
db.execute('SET threads=8')
raw='D:/Brazil_Acute_Inflammation_Climate_2003_2023/hybrid_work_v1/private/health_partitions/rd/*/*/*.parquet'
xw='D:/Brazil_Acute_Inflammation_Climate_2003_2023/hybrid_work_v1/reference_step03/municipality_sih6_ibge7_health_region_2023_v2023_1.csv'
db.execute(f"""CREATE TABLE carry AS SELECT x.health_region_2023::VARCHAR hrid,
 try_strptime(r.DT_INTER,'%Y%m%d')::DATE d1,try_strptime(r.DT_SAIDA,'%Y%m%d')::DATE d2,
 count(*) admissions FROM read_parquet('{raw}') r
 LEFT JOIN read_csv_auto('{xw}') x ON r.MUNIC_RES=x.municipality_sih6::VARCHAR
 WHERE substr(r.DIAG_PRINC,1,3)='N20' AND r.CAR_INT='02'
 AND r.DT_INTER<'20180101' AND r.DT_SAIDA>='20180101'
 GROUP BY ALL""")
report=db.sql("SELECT sum(admissions) carry_in_admissions,sum(admissions*(date_diff('day',DATE '2018-01-01',least(d2,DATE '2022-12-31'))+1)) carry_in_inclusive_days,min(d1) earliest_admission,max(d2) latest_discharge FROM carry").df()
report.to_csv(out/'occupancy_boundary_audit.csv',index=False,encoding='utf-8-sig')
print(report.to_dict('records'),flush=True)
db.execute("CREATE TABLE intervals AS SELECT hrid,d1,d2,admissions FROM read_parquet('"+(work/'deidentified_interval_counts.parquet').as_posix()+"') UNION ALL SELECT * FROM carry")
db.execute("""CREATE TABLE grid AS SELECT DISTINCT hrid,unnest(generate_series(DATE '2018-01-01',DATE '2022-12-31',INTERVAL 1 DAY))::DATE date FROM intervals""")
for label,exit_delta in [('inclusive','1'),('exclusive','0')]:
    db.execute(f"""CREATE OR REPLACE TABLE changes AS
     SELECT hrid,greatest(d1,DATE '2018-01-01') AS date,sum(admissions) change
     FROM intervals WHERE d2>=d1 {'AND d2>d1' if label=='exclusive' else ''} GROUP BY 1,2
     UNION ALL SELECT hrid,d2+{exit_delta}::INTEGER AS date,-sum(admissions) change
     FROM intervals WHERE d2>=d1 {'AND d2>d1' if label=='exclusive' else ''} GROUP BY 1,2""")
    db.sql("""WITH c AS(SELECT hrid,date,sum(change) change FROM changes GROUP BY 1,2)
     SELECT g.hrid,g.date,sum(coalesce(c.change,0)) OVER(PARTITION BY g.hrid ORDER BY g.date) occupied
     FROM grid g LEFT JOIN c USING(hrid,date)""").write_parquet(str(work/f'occupied_boundary_{label}.parquet'))
print('BOUNDARY_COMPLETE',flush=True)
