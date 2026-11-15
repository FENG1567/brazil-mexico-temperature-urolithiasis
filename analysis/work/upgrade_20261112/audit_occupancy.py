"""Audit interval reconstruction using aggregate durations only."""
import sys
from pathlib import Path
import json
sys.path.insert(0,'F:/analysis/climate_biliary_20261102/work/pylibs')
import duckdb
base=Path(__file__).resolve().parents[2]
out=base/'outputs/upgrade_20261112'
work=base/'work/upgrade_20261112'
db=duckdb.connect()
db.execute('SET threads=8')
raw='D:/Brazil_Acute_Inflammation_Climate_2003_2023/hybrid_work_v1/private/health_partitions/rd/*/*/*.parquet'
xw='D:/Brazil_Acute_Inflammation_Climate_2003_2023/hybrid_work_v1/reference_step03/municipality_sih6_ibge7_health_region_2023_v2023_1.csv'
db.execute(f"""CREATE TABLE intervals AS SELECT x.health_region_2023::VARCHAR hrid,
 try_strptime(replace(r.DT_INTER,'-',''),'%Y%m%d')::DATE d1,
 try_strptime(replace(r.DT_SAIDA,'-',''),'%Y%m%d')::DATE d2,
 r.DIAS_PERM beddays,count(*) admissions
 FROM read_parquet('{raw}') r
 LEFT JOIN read_csv_auto('{xw}') x ON r.MUNIC_RES=x.municipality_sih6::VARCHAR
 WHERE substr(r.DIAG_PRINC,1,3)='N20' AND r.CAR_INT='02'
 AND substr(replace(r.DT_INTER,'-',''),1,4) BETWEEN '2018' AND '2022'
 GROUP BY ALL""")
report=db.sql("""SELECT sum(admissions) admissions,
 sum(admissions) FILTER(WHERE d1 IS NULL OR d2 IS NULL) missing_date,
 sum(admissions) FILTER(WHERE d2<d1) reversed_interval,
 sum(admissions) FILTER(WHERE hrid IS NULL) unmapped,
 sum(admissions) FILTER(WHERE d1=d2) same_day,
 sum(admissions) FILTER(WHERE d2>DATE '2022-12-31') discharged_after_window,
 sum(admissions*beddays) recorded_beddays,
 sum(admissions*date_diff('day',d1,d2)) elapsed_days,
 sum(admissions*(date_diff('day',d1,d2)+1)) inclusive_days,
 sum(admissions) FILTER(WHERE beddays<>date_diff('day',d1,d2)) duration_disagreement
 FROM intervals""").df()
report=report.fillna(0)
report.to_csv(out/'occupancy_interval_audit.csv',index=False,encoding='utf-8-sig')
print(report.to_dict('records'),flush=True)
db.sql('SELECT * FROM intervals').write_parquet(str(work/'deidentified_interval_counts.parquet'))
# Discharge-exclusive occupancy is an explicit alternative to the old inclusive convention.
db.execute("""CREATE TABLE delta AS SELECT hrid,d1 AS date,sum(admissions) change
 FROM intervals WHERE d2>d1 GROUP BY 1,2
 UNION ALL SELECT hrid,d2 AS date,-sum(admissions) change
 FROM intervals WHERE d2>d1 GROUP BY 1,2""")
db.sql("""WITH changes AS (SELECT hrid,date,sum(change) change FROM delta GROUP BY 1,2),
 grid AS (SELECT DISTINCT hrid,unnest(generate_series(DATE '2018-01-01',DATE '2022-12-31',INTERVAL 1 DAY))::DATE date FROM intervals)
 SELECT g.hrid,g.date,sum(coalesce(c.change,0)) OVER(PARTITION BY g.hrid ORDER BY g.date) occupied_exclusive
 FROM grid g LEFT JOIN changes c USING(hrid,date)""").write_parquet(str(work/'occupied_discharge_exclusive.parquet'))
print('OCCUPANCY_AUDIT_COMPLETE',flush=True)
