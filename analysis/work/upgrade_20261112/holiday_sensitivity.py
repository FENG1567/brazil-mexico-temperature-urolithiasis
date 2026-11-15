"""Six calendar-control sensitivity models specified before these fits, 2018–2022.

The calendar is generated with python-holidays' PUBLIC category and
observed=True, checked against recurring date rules, and supplemented with
nationwide federal election dates in Brazil and Mexico as operational controls.
These are additional sensitivities after the earlier main results were known.
All outcome inputs are existing de-identified state/day or region/day frames.
"""

import os

for variable in ("OMP_NUM_THREADS", "NUMBA_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[variable] = "8"

import hashlib
import importlib.metadata
import importlib.util
import json
import sys
import time
from datetime import date, timedelta
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
WORK = BASE / "work/upgrade_20261112"
OUT = BASE / "outputs/upgrade_20261112"
CALENDAR_LIBS = BASE / "work/calendarlibs"
SOURCE_WORK = Path("F:/analysis/climate_biliary_20261102/work")
CLINICAL_PATH = WORK / "clinical_code_day_aggregate.parquet"
SP_LINK_PATH = OUT / "sp_linked_daily_outcomes.csv"
PANEL_PATH = SOURCE_WORK / "panel_composite_humidity.parquet"
YEARS = list(range(2018, 2023))
GROUPS = ["neither", "heat_only", "smoke_only", "both"]
COUNTRY_NAMES = {"BR": "brazil", "MX": "mexico"}
STATE_STRICT_TOTALS = {"BR": 321_201, "MX": 32_417}
BR_URGENT_ELECTIVE = {"broad_urgent": 315_504, "broad_elective": 115_870,
                      "strict_urgent": 226_402, "strict_elective": 94_799}
SP_EVENT_MINIMUMS = {"any_drainage": 47_918, "any_removal": 35_976}

sys.path.insert(0, str(CALENDAR_LIBS))
sys.path.insert(0, "F:/analysis/climate_biliary_20261102/work/pylibs")

import duckdb
import holidays
import numpy as np
import pandas as pd
import pyfixest as pf
from holidays.constants import PUBLIC
from scipy.stats import norm, t as student_t


OUT.mkdir(parents=True, exist_ok=True)


OFFICIAL_SOURCES = [
    {
        "country": "BR",
        "title": "Lei 662/1949, updated by Lei 10.607/2002",
        "url": "https://www2.camara.leg.br/legin/fed/lei/1940-1949/lei-662-6-abril-1949-347136-normaatualizada-pl.html",
        "supports": "National fixed holidays, including Tiradentes and Finados in the current wording.",
    },
    {
        "country": "BR",
        "title": "Lei 6.802/1980",
        "url": "https://legis.senado.gov.br/norma/548678/publicacao/15645366",
        "supports": "12 October as a national holiday.",
    },
    {
        "country": "BR",
        "title": "Lei 9.093/1995",
        "url": "https://www.planalto.gov.br/ccivil_03/leis/l9093.htm",
        "supports": "Framework for civil and religious holidays; Good Friday is recognized in annual federal holiday calendars.",
    },
    {
        "country": "BR",
        "title": "TCU Portaria 28/2018",
        "url": "https://pesquisa.apps.tcu.gov.br/doc/norma/Portaria/TCU/28/2018/",
        "supports": "2018 federal holiday calendar and Good Friday date.",
    },
    {
        "country": "BR",
        "title": "Ministry of Planning Portaria 444/2018 for 2019",
        "url": "https://www.gov.br/economia/pt-br/assuntos/planejamento/orcamento/orcamentos-anuais/2018/legislacao/programacao/portaria-mp-no-444-de-27-de-dezembro-de-2018.pdf",
        "supports": "2019 federal holiday calendar and Good Friday date.",
    },
    {
        "country": "BR",
        "title": "Federal District 2020 official holiday calendar",
        "url": "https://www.economia.df.gov.br/calendario-de-feriados-e-pontos-facultativos-de-2020",
        "supports": "2020 national holiday dates, including Good Friday; local and optional days are excluded from this calendar.",
    },
    {
        "country": "BR",
        "title": "TCU Portaria 14/2021",
        "url": "https://pesquisa.apps.tcu.gov.br/doc/norma/Portaria/TCU/14/2021/",
        "supports": "2021 federal holiday calendar and Good Friday date.",
    },
    {
        "country": "BR",
        "title": "Federal Government 2022 holiday calendar",
        "url": "https://www.gov.br/pt-br/noticias/financas-impostos-e-gestao-publica/2021/12/governo-federal-divulga-calendario-2022-de-feriados-e-pontos-facultativos",
        "supports": "2022 federal holiday calendar and Good Friday date.",
    },
    {
        "country": "BR",
        "title": "Lei 10.607/2002, revocation of Lei 1.266/1950",
        "url": "https://www2.camara.leg.br/legin/fed/lei/2002/lei-10607-19-dezembro-2002-473572-publicacaooriginal-1-pl.html",
        "supports": "Revokes the earlier statute declaring nationwide election days holidays. Brazilian election dates are included separately as operational calendar controls, without assigning a statutory holiday classification.",
    },
    {
        "country": "MX",
        "title": "Ley Federal del Trabajo, Article 74, last amendment through 2021",
        "url": "https://portalhcd.diputados.gob.mx/LeyesBiblio/PortalWeb/Leyes/Vigentes/PDF/125_230421.pdf",
        "supports": "Mandatory rest days, Monday-observed February/March/November dates, six-year executive transition day, and ordinary election day.",
    },
    {
        "country": "MX",
        "title": "INE 2018 federal election calendar",
        "url": "https://www.ine.mx/voto-y-elecciones/elecciones-2018/eleccion-federal/",
        "supports": "Federal election day, 1 July 2018; added under LFT Article 74(IX).",
    },
    {
        "country": "MX",
        "title": "INE 2021 federal election calendar",
        "url": "https://www.ine.mx/voto-y-elecciones/elecciones-2021/eleccion-federal-2021/",
        "supports": "Federal election day, 6 June 2021; added under LFT Article 74(IX).",
    },
]

SPECIFICATIONS = [
    {
        "model": "brazil_strict_n20_holiday_lag07",
        "sample": "Brazil 27 states; strict N20 all admissions; 2018-2022",
        "formula": "strict_N20 ~ t0 + t1 + t2 + t3 + t4 + t5 + t6 + t7 + holiday | state_year_month_day_of_week",
        "cluster": "two-way CRV1 by state and date; t inference with df=min(G_state,G_date)-1",
    },
    {
        "model": "mexico_strict_n20_holiday_lag07",
        "sample": "Mexico 32 states; strict N20 all admissions; ERA5-Land temperature; 2018-2022",
        "formula": "strict_N20 ~ t0 + t1 + t2 + t3 + t4 + t5 + t6 + t7 + holiday | state_year_month_day_of_week",
        "cluster": "two-way CRV1 by state and date; t inference with df=min(G_state,G_date)-1",
    },
    {
        "model": "sp_any_drainage_holiday_lag07",
        "sample": "Brazil health-region days; any drainage outcome; 2018-2022",
        "formula": "any_drainage ~ t0 + t1 + t2 + t3 + t4 + t5 + t6 + t7 + f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + holiday | hrid_year_month_day_of_week",
        "cluster": "CRV1 by hrid",
    },
    {
        "model": "sp_any_removal_holiday_lag07",
        "sample": "Brazil health-region days; any removal outcome; 2018-2022",
        "formula": "any_removal ~ t0 + t1 + t2 + t3 + t4 + t5 + t6 + t7 + f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + holiday | hrid_year_month_day_of_week",
        "cluster": "CRV1 by hrid",
    },
    {
        "model": "broad_urgent_elective_holiday_joint",
        "sample": "Brazil region days stacked as N20-N23 urgent/elective counts; 2018-2022",
        "formula": "count ~ C(outcome)*C(group) + C(outcome):holiday | stratum^outcome",
        "cluster": "CRV1 by hrid",
    },
    {
        "model": "strict_urgent_elective_holiday_joint",
        "sample": "Brazil region days stacked as strict N20 urgent/elective counts; 2018-2022",
        "formula": "count ~ C(outcome)*C(group) + C(outcome):holiday | stratum^outcome",
        "cluster": "CRV1 by hrid",
    },
]


def sql_path(path):
    return Path(path).as_posix().replace("'", "''")


def import_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load existing analysis helper: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def first_weekday(year, month, weekday=0):
    day = date(year, month, 1)
    return day + timedelta(days=(weekday - day.weekday()) % 7)


def nth_weekday(year, month, n, weekday=0):
    return first_weekday(year, month, weekday) + timedelta(days=7 * (n - 1))


def official_expected_dates():
    """Encode the official national holiday rules independently of the package."""
    brazil = set()
    for year in YEARS:
        brazil.update({
            date(year, 1, 1), date(year, 4, 21), date(year, 5, 1),
            date(year, 9, 7), date(year, 10, 12), date(year, 11, 2),
            date(year, 11, 15), date(year, 12, 25),
        })
        # Federal holiday calendars list Paixão de Cristo each year.
        import dateutil.easter
        brazil.add(dateutil.easter.easter(year) - timedelta(days=2))

    mexico = set()
    for year in YEARS:
        mexico.update({
            date(year, 1, 1),
            nth_weekday(year, 2, 1),
            nth_weekday(year, 3, 3),
            date(year, 5, 1),
            date(year, 9, 16),
            nth_weekday(year, 11, 3),
            date(year, 12, 25),
        })
        if year == 2018:
            mexico.add(date(year, 12, 1))
    return {"BR": brazil, "MX": mexico}


def observed_role(country, holiday_date, name):
    if country == "MX":
        historical = {
            "Día de la Constitución": (2, 5),
            "Natalicio de Benito Juárez": (3, 21),
            "Día de la Revolución": (11, 20),
        }
        for title, (month, day) in historical.items():
            if title in name and (holiday_date.month, holiday_date.day) != (month, day):
                return "statutory weekday-observed date"
    return "statutory calendar date"


def build_holiday_calendar():
    package_version = importlib.metadata.version("holidays")
    if package_version != "0.108":
        raise ValueError(f"Expected pinned holidays version 0.108, found {package_version}")
    expected = official_expected_dates()
    records = []
    package_dates = {}
    for country, language in (("BR", "pt_BR"), ("MX", "es")):
        calendar = holidays.country_holidays(
            country,
            years=YEARS,
            observed=True,
            categories=PUBLIC,
            language=language,
        )
        actual_dates = {day for day in calendar.keys() if day.year in YEARS}
        package_dates[country] = actual_dates
        if actual_dates != expected[country]:
            missing = sorted(expected[country] - actual_dates)
            extra = sorted(actual_dates - expected[country])
            raise ValueError(f"{country} library calendar differs from official date rules; missing={missing}; extra={extra}")
        for holiday_date, name in sorted(calendar.items()):
            if holiday_date.year not in YEARS:
                continue
            records.append({
                "country": country,
                "date": holiday_date.isoformat(),
                "holiday_name": str(name),
                "calendar_source": "python-holidays PUBLIC",
                "library_version": package_version,
                "observed_argument": True,
                "date_role": observed_role(country, holiday_date, str(name)),
            })

    # Add nationwide federal elections as operational calendar controls in both
    # countries. This definition avoids inferring hospital schedules from the
    # legal classification of a Brazilian election date.
    election_days = [
        ("BR", "2018-10-07", "Eleições gerais 2018 primeiro turno"),
        ("BR", "2018-10-28", "Eleições gerais 2018 segundo turno"),
        ("BR", "2022-10-02", "Eleições gerais 2022 primeiro turno"),
        ("BR", "2022-10-30", "Eleições gerais 2022 segundo turno"),
        ("MX", "2018-07-01", "Jornada electoral federal 2018"),
        ("MX", "2021-06-06", "Jornada electoral federal 2021"),
    ]
    for country, holiday_date, name in election_days:
        records.append({
            "country": country,
            "date": holiday_date,
            "holiday_name": name,
            "calendar_source": "Official TSE federal election calendar" if country == "BR" else "LFT Article 74(IX) + official INE federal calendar",
            "library_version": package_version,
            "observed_argument": True,
            "date_role": "nationwide ordinary federal election day",
        })
    calendar_frame = pd.DataFrame(records).sort_values(["country", "date", "holiday_name"]).reset_index(drop=True)
    if calendar_frame.duplicated(["country", "date"]).any():
        raise ValueError("National holiday calendar has duplicate country/date entries")

    core_counts = {country: len(package_dates[country]) for country in ("BR", "MX")}
    all_counts = calendar_frame.groupby("country").size().astype(int).to_dict()
    provenance = {
        "calendar_window": "2018-01-01 through 2022-12-31",
        "python_holidays_version": package_version,
        "pypi_project": "https://pypi.org/project/holidays/0.108/",
        "package_location": str(CALENDAR_LIBS),
        "category": "PUBLIC",
        "observed": True,
        "language": {"BR": "pt_BR", "MX": "es"},
        "source_method": "Generate recurring PUBLIC dates with python-holidays, observed=True, verify the recurring rules, and add officially dated nationwide federal elections in both countries as operational calendar controls.",
        "official_rule_checks": {
            "BR": "Lei 662/1949 as amended by 10.607/2002, Lei 6.802/1980; annual federal holiday calendars confirm Good Friday dates. Fixed statutory dates remain on their legal dates; no national substitute date is added by the package.",
            "MX": "Ley Federal del Trabajo Article 74(I-VIII) as in force through 2022; first Monday in February, third Monday in March, third Monday in November, and 1 December 2018 executive transition date are matched exactly.",
        },
        "supplemental_dates": [
            {"country": country, "date": day, "name": name,
             "basis": "Nationwide federal election date confirmed by TSE; operational calendar control." if country == "BR" else "LFT Article 74(IX), nationwide federal ordinary election; date confirmed by INE official calendar."}
            for country, day, name in election_days
        ],
        "date_counts": {"python_holidays_core": core_counts, "after_national_election_supplement": all_counts},
        "included": "Recurring PUBLIC calendar dates, observed Monday dates in Mexico, and nationwide federal election days in both countries.",
        "excluded": "State/local holidays, municipal religious holidays beyond the recurring country calendar, administrative optional days, Carnival/bridge days outside the PUBLIC category, and state/local-only elections.",
        "known_coverage_limit": "This indicator does not represent all local schedules, hospital-specific staffing, bridge days, or every state/local election-related closure.",
        "official_sources": OFFICIAL_SOURCES + [
            {"title": "TSE chronology of elections, 2018 dates", "url": "https://www.tse.jus.br/eleicoes/historia/cronologia-das-eleicoes"},
            {"title": "TSE official 2022 election calendar", "url": "https://www.tse.jus.br/comunicacao/noticias/2022/Janeiro/confira-as-principais-datas-do-calendario-eleitoral-de-2022"},
        ],
    }
    calendar_frame.to_csv(OUT / "holiday_calendar_dates.csv", index=False, encoding="utf-8-sig")
    (OUT / "holiday_calendar_provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8"
    )
    return calendar_frame, provenance


def add_holiday_flag(frame, country_column="country", date_column="date", calendar_frame=None):
    keys = set(zip(calendar_frame.country, pd.to_datetime(calendar_frame.date).dt.normalize()))
    return np.fromiter(
        ((country, pd.Timestamp(day).normalize()) in keys
         for country, day in zip(frame[country_column], frame[date_column])),
        dtype=np.int8,
        count=len(frame),
    )


def estimate(weights, names, coef, covariance, df=None, inference="normal"):
    vector = np.zeros(len(names), dtype=float)
    for term, weight in weights.items():
        if term not in names:
            raise ValueError(f"Required contrast term missing: {term}; available={names}")
        vector[names.index(term)] = weight
    beta = float(vector @ coef)
    variance = float(vector @ covariance @ vector)
    if variance < 0 and abs(variance) < 1e-12:
        variance = 0.0
    if variance < 0:
        raise ValueError(f"Negative contrast variance: {variance}")
    se = float(np.sqrt(variance))
    if inference == "t":
        if df is None or df <= 0:
            raise ValueError(f"Invalid t degrees of freedom: {df}")
        critical = float(student_t.ppf(0.975, df))
        p = float(2 * student_t.sf(abs(beta / se), df)) if se else 0.0
    else:
        critical = 1.96
        p = float(2 * norm.sf(abs(beta / se))) if se else 0.0
    return {
        "beta": beta,
        "variance": variance,
        "se": se,
        "RR": float(np.exp(beta)),
        "CI_low": float(np.exp(beta - critical * se)),
        "CI_high": float(np.exp(beta + critical * se)),
        "p": p,
        "inference": inference,
        "df": df if inference == "t" else None,
        "weights": vector.tolist(),
        "contrast_terms": [term for term, weight in weights.items() if weight != 0],
    }


def fit_and_save(model_key, formula, frame, outcome, clusters, cluster_fields,
                 contrasts, metadata):
    started = time.time()
    print(f"START {model_key} rows={len(frame):,} formula={formula}", flush=True)
    model = pf.fepois(formula, data=frame, vcov={"CRV1": clusters})
    coef_series = model.coef()
    names = [str(name) for name in coef_series.index]
    coef = coef_series.to_numpy(dtype=float)
    covariance = np.asarray(model._vcov, dtype=float)
    if covariance.shape != (len(names), len(names)):
        raise ValueError(f"Unexpected covariance shape for {model_key}: {covariance.shape}")

    results = []
    coefficients = []
    for position, term in enumerate(names):
        term_estimate = estimate(
            {term: 1.0}, names, coef, covariance,
            df=metadata.get("inference_df"), inference=metadata.get("inference", "normal"),
        )
        coefficients.append({
            "model": model_key,
            "term": term,
            "beta": term_estimate["beta"],
            "se": term_estimate["se"],
            "variance": float(covariance[position, position]),
            "RR": term_estimate["RR"],
            "CI_low": term_estimate["CI_low"],
            "CI_high": term_estimate["CI_high"],
            "p": term_estimate["p"],
            "sign": "positive" if term_estimate["beta"] > 0 else "negative" if term_estimate["beta"] < 0 else "zero",
            "inference": term_estimate["inference"],
            "df": term_estimate["df"],
            "formula": formula,
        })

    for contrast_name, weights in contrasts.items():
        result = estimate(
            weights,
            names,
            coef,
            covariance,
            df=metadata.get("inference_df"),
            inference=metadata.get("inference", "normal"),
        )
        results.append({"model": model_key, "contrast": contrast_name, **result})

    observed_outcome = float(frame[outcome].sum()) if outcome in frame.columns else float(frame["count"].sum())
    summary = {
        "model": model_key,
        "formula": formula,
        "cluster": f"CRV1 by {clusters}",
        "inference": metadata.get("inference", "normal"),
        "inference_df": metadata.get("inference_df"),
        "input_rows": int(len(frame)),
        "model_rows": int(model._N),
        "rows_dropped_by_fit": int(len(frame) - int(model._N)),
        "input_outcome_events": observed_outcome,
        "model_outcome_events": float(model._data[outcome].sum()) if outcome in frame.columns else None,
        "cluster_counts": {field: int(frame[field].nunique()) for field in cluster_fields},
        "strata_input": int(frame.stratum.nunique()) if "stratum" in frame.columns else None,
        "holiday_input_rows": int(frame.holiday.sum()),
        "holiday_model_rows": None,
        "seconds": round(time.time() - started, 2),
        **metadata.get("summary_fields", {}),
    }
    details = {
        **summary,
        "coefficient_names": names,
        "coefficients": coef.tolist(),
        "covariance": covariance.tolist(),
        "contrasts": {row["contrast"]: row for row in results},
    }
    print(
        f"FIT_OK {model_key} N={summary['model_rows']:,} "
        f"events={summary['model_outcome_events']} terms={names}",
        flush=True,
    )
    return coefficients, results, summary, details


def state_level_frame(mexico_module):
    mexico, mexico_era5, _ = mexico_module.load_mexico()
    brazil = mexico_module.load_brazil()
    mexico_state = mexico[["state", "date", "strict_N20"]].merge(
        mexico_era5[["state", "date", "temperature"]],
        on=["state", "date"], how="left", validate="one_to_one",
    )
    mexico_state["country"] = "MX"
    mexico_state["temperature_source"] = "ERA5-Land area-weighted state daily"
    brazil_state = brazil[["state", "date", "strict_N20", "tmean"]].rename(
        columns={"tmean": "temperature"}
    )
    brazil_state["country"] = "BR"
    brazil_state["temperature_source"] = "validated Brazil state population-weighted panel temperature"
    frame = pd.concat([brazil_state, mexico_state], ignore_index=True)
    frame["state"] = frame.state.astype(str).str.zfill(2)
    frame["date"] = pd.to_datetime(frame.date).dt.normalize()
    frame["country_state"] = frame.country + "_" + frame.state
    frame["date_id"] = frame.date.dt.strftime("%Y%m%d")
    frame["year_month"] = frame.date.dt.strftime("%Y%m")
    frame["day_of_week"] = frame.date.dt.day_of_week.astype(int)
    frame["stratum"] = frame.country_state + "_" + frame.year_month + "_" + frame.day_of_week.astype(str)
    frame = frame.sort_values(["country_state", "date"]).reset_index(drop=True)
    if len(frame) != 59 * 1826 or frame.country_state.nunique() != 59 or frame.duplicated(["country_state", "date"]).any():
        raise ValueError(f"Expected complete 59-state combined calendar; got {len(frame)} rows")
    totals = frame.groupby("country", observed=True).strict_N20.sum().astype(int).to_dict()
    if totals != STATE_STRICT_TOTALS:
        raise ValueError(f"Strict state-level event totals changed: {totals}")
    return frame


def fit_state_strict(frame, country, calendar):
    country_name = COUNTRY_NAMES[country]
    frame = frame.copy()
    frame["holiday"] = add_holiday_flag(frame, calendar_frame=calendar)
    for lag in range(8):
        frame[f"t{lag}"] = frame.groupby("country_state", sort=False).temperature.shift(lag)
    lag_terms = [f"t{lag}" for lag in range(8)]
    model_frame = frame.dropna(subset=lag_terms).copy()
    expected_states = {"BR": 27, "MX": 32}[country]
    if frame.country_state.nunique() != expected_states or len(frame) != expected_states * 1826:
        raise ValueError(f"Unexpected {country} state calendar dimensions")
    if len(model_frame) != expected_states * (1826 - 7):
        raise ValueError(f"Unexpected state-lag sample size: {len(model_frame)}")
    state_clusters = model_frame.country_state.nunique()
    date_clusters = model_frame.date_id.nunique()
    inference_df = min(state_clusters, date_clusters) - 1
    formula = "strict_N20 ~ " + " + ".join(lag_terms + ["holiday"]) + " | stratum"
    terms = {term: 1.0 for term in lag_terms}
    contrasts = {
        "cumulative_temperature_lag_0_7_per_C": terms,
        "national_holiday_rate_ratio": {"holiday": 1.0},
    }
    return fit_and_save(
        f"{country_name}_strict_n20_holiday_lag07",
        formula,
        model_frame,
        "strict_N20",
        "country_state + date_id",
        ["country_state", "date_id"],
        contrasts,
        {
            "inference": "t",
            "inference_df": inference_df,
            "summary_fields": {
                "study_states": int(frame.country_state.nunique()),
                "country": country,
                "strict_N20_events": int(frame.strict_N20.sum()),
                "input_days_after_lag_warmup": int(len(model_frame)),
                "national_holiday_days_input": int(model_frame.holiday.sum()),
                "temperature_sources": sorted(frame.temperature_source.unique().tolist()),
            },
        },
    )


def load_sp_frame(sp_module, calendar):
    # prepare_panel() is a validated helper; suppress its unrelated audit write.
    sp_module.save_json = lambda *_args, **_kwargs: None
    frame, audit = sp_module.prepare_panel()
    frame = frame.copy()
    frame["date"] = pd.to_datetime(frame.date).dt.normalize()
    frame["holiday"] = add_holiday_flag(
        frame.assign(country="BR"), country_column="country", calendar_frame=calendar
    )
    for lag in range(8):
        frame[f"t{lag}"] = frame.groupby("hrid", sort=False).tmean.shift(lag)
        frame[f"f{lag}"] = frame.groupby("hrid", sort=False).fire10.shift(lag)
    lag_terms = [f"t{lag}" for lag in range(8)] + [f"f{lag}" for lag in range(8)]
    model_frame = frame.dropna(subset=lag_terms).copy()
    if len(frame) != 821_700 or len(model_frame) != 818_550:
        raise ValueError(f"Unexpected SP continuous-lag sample: full={len(frame)}, lagged={len(model_frame)}")
    for outcome, minimum in SP_EVENT_MINIMUMS.items():
        if int(model_frame[outcome].sum()) < minimum:
            raise ValueError(f"Validated SP outcome total changed for {outcome}")
    return frame, model_frame, lag_terms, audit


def fit_sp_sensitivities(full_frame, model_frame, lag_terms, audit, on_fit=None):
    coefficients, results, summaries, models = [], [], [], {}
    for outcome in ("any_drainage", "any_removal"):
        formula = outcome + " ~ " + " + ".join(lag_terms + ["holiday"]) + " | stratum"
        temperature_weights = {f"t{lag}": 1.0 for lag in range(8)}
        fire_weights = {f"f{lag}": 1.0 for lag in range(8)}
        contrasts = {
            "cumulative_temperature_lag_0_7_per_C": temperature_weights,
            "cumulative_fire_PM25_lag_0_7_per_10_units": fire_weights,
            "national_holiday_rate_ratio": {"holiday": 1.0},
        }
        coef_rows, contrast_rows, summary, model = fit_and_save(
            f"sp_{outcome}_holiday_lag07",
            formula,
            model_frame,
            outcome,
            "hrid",
            ["hrid"],
            contrasts,
            {
                "inference": "normal",
                "summary_fields": {
                    "analysis": "preplanned additional SP-linked sensitivity",
                    "sp_source_sha256": audit["SP_input_sha256"],
                    "source_full_window_events": int(full_frame[outcome].sum()),
                    "fire_scale": "fire_pm25 / 10",
                    "input_days_after_lag_warmup": int(len(model_frame)),
                },
            },
        )
        coefficients.extend(coef_rows)
        results.extend(contrast_rows)
        summaries.append(summary)
        models[summary["model"]] = model
        if on_fit is not None:
            on_fit(coef_rows, contrast_rows, summary, model)
    return coefficients, results, summaries, models


def brazil_urgent_elective_panel(sp_frame):
    db = duckdb.connect()
    db.execute("SET threads=8")
    path = sql_path(CLINICAL_PATH)
    daily = db.sql(
        f"""
        SELECT CAST(hrid AS VARCHAR) AS hrid, CAST(date AS DATE) AS date,
               SUM(CASE WHEN diagnosis IN ('N20','N21','N22','N23') AND CAR_INT='02'
                        THEN admissions ELSE 0 END) AS broad_urgent,
               SUM(CASE WHEN diagnosis IN ('N20','N21','N22','N23') AND CAR_INT='01'
                        THEN admissions ELSE 0 END) AS broad_elective,
               SUM(CASE WHEN diagnosis='N20' AND CAR_INT='02'
                        THEN admissions ELSE 0 END) AS strict_urgent,
               SUM(CASE WHEN diagnosis='N20' AND CAR_INT='01'
                        THEN admissions ELSE 0 END) AS strict_elective
        FROM read_parquet('{path}')
        WHERE IDENT='1' AND date >= DATE '2018-01-01' AND date < DATE '2023-01-01'
        GROUP BY 1,2
        """
    ).df()
    db.close()
    daily["hrid"] = daily.hrid.astype(str)
    daily["date"] = pd.to_datetime(daily.date).dt.normalize()
    frame = sp_frame.loc[sp_frame.date.between("2018-01-01", "2022-12-31")].copy()
    for column in ("broad_urgent", "broad_elective", "strict_urgent", "strict_elective"):
        frame = frame.drop(columns=column, errors="ignore")
    frame = frame.merge(daily, on=["hrid", "date"], how="left", validate="one_to_one")
    endpoints = list(BR_URGENT_ELECTIVE)
    frame[endpoints] = frame[endpoints].fillna(0).astype("int64")
    totals = frame[endpoints].sum().astype(int).to_dict()
    if totals != BR_URGENT_ELECTIVE:
        raise ValueError(f"Brazil urgent/elective source counts differ: {totals}")
    frame["group"] = np.select(
        [
            frame.heat.eq(0) & frame.smoke.eq(0),
            frame.heat.eq(1) & frame.smoke.eq(0),
            frame.heat.eq(0) & frame.smoke.eq(1),
            frame.heat.eq(1) & frame.smoke.eq(1),
        ],
        GROUPS,
        default="other",
    )
    if frame.group.eq("other").any():
        raise ValueError("Brazil urgent/elective panel has an out-of-scope exposure category")
    frame["group"] = pd.Categorical(frame.group, categories=GROUPS)
    frame["country"] = "BR"
    frame["holiday"] = add_holiday_flag(frame, calendar_frame=holiday_calendar_global[0])
    return frame


holiday_calendar_global = [None]


def fit_joint_urgent_elective(frame, phenotype):
    if phenotype == "broad":
        endpoints = ["broad_elective", "broad_urgent"]
        key = "broad_urgent_elective_holiday_joint"
        event_columns = {"elective": "broad_elective", "urgent": "broad_urgent"}
    else:
        endpoints = ["strict_elective", "strict_urgent"]
        key = "strict_urgent_elective_holiday_joint"
        event_columns = {"elective": "strict_elective", "urgent": "strict_urgent"}

    long = frame[["hrid", "date", "stratum", "group", "holiday", *endpoints]].melt(
        id_vars=["hrid", "date", "stratum", "group", "holiday"],
        value_vars=endpoints,
        var_name="outcome_source",
        value_name="count",
    )
    source_to_outcome = {source: label for label, source in event_columns.items()}
    long["outcome"] = pd.Categorical(
        long.outcome_source.map(source_to_outcome), categories=["elective", "urgent"]
    )
    long["endpoint_stratum"] = long.stratum + "_" + long.outcome.astype(str)
    formula = "count ~ C(outcome)*C(group) + C(outcome):holiday | stratum^outcome"
    contrasts = {
        f"{group}_urgent_effect_over_elective_effect": {
            f"C(outcome)[T.urgent]:C(group)[T.{group}]": 1.0
        }
        for group in GROUPS[1:]
    }
    coefficients, results, summary, details = fit_and_save(
        key,
        formula,
        long,
        "count",
        "hrid",
        ["hrid"],
        contrasts,
        {
            "inference": "normal",
            "summary_fields": {
                "phenotype": phenotype,
                "stacked_endpoints": ["elective", "urgent"],
                "stratum_definition": "health region × year-month × day-of-week, separately interacted with endpoint",
                "holiday_coefficient_definition": "separate C(outcome):holiday coefficients for elective and urgent",
                "input_days_per_endpoint": int(len(frame)),
                "input_events_by_endpoint": {
                    label: int(long.loc[long.outcome.eq(label), "count"].sum())
                    for label in ("elective", "urgent")
                },
            },
        },
    )

    names = details["coefficient_names"]
    coef = np.asarray(details["coefficients"], dtype=float)
    covariance = np.asarray(details["covariance"], dtype=float)
    holiday_terms = [name for name in names if name.endswith(":holiday")]
    elective_holiday_term = next((name for name in holiday_terms if "elective" in name), None)
    urgent_holiday_term = next((name for name in holiday_terms if "urgent" in name), None)
    if elective_holiday_term is None or urgent_holiday_term is None:
        raise ValueError(f"Expected endpoint-specific holiday terms, found {holiday_terms}")
    endpoint_holiday = {
        "elective_holiday_rate_ratio": {elective_holiday_term: 1.0},
        "urgent_holiday_rate_ratio": {urgent_holiday_term: 1.0},
        "urgent_vs_elective_holiday_rate_ratio": {
            urgent_holiday_term: 1.0, elective_holiday_term: -1.0
        },
    }
    for name, weights in endpoint_holiday.items():
        estimate_row = estimate(weights, names, coef, covariance)
        results.append({"model": key, "contrast": name, **estimate_row})
    details["contrasts"].update({row["contrast"]: row for row in results if row["model"] == key})
    for result in results:
        if result["model"] == key and result["contrast"] in endpoint_holiday:
            coefficients.append({
                "model": key,
                "term": result["contrast"],
                "beta": result["beta"],
                "se": result["se"],
                "variance": result["variance"],
                "RR": result["RR"],
                "CI_low": result["CI_low"],
                "CI_high": result["CI_high"],
                "p": result["p"],
                "sign": "positive" if result["beta"] > 0 else "negative" if result["beta"] < 0 else "zero",
                "inference": "normal",
                "df": None,
                "formula": formula,
            })
    return coefficients, results, summary, details


def make_report(calendar_frame, state_results, sp_results, joint_results,
                summaries, provenance):
    lookup = {(row["model"], row["contrast"]): row for row in state_results + sp_results + joint_results}
    summary_lookup = {row["model"]: row for row in summaries}

    def fmt(model, contrast):
        row = lookup[(model, contrast)]
        return f"{row['RR']:.4f} ({row['CI_low']:.4f}–{row['CI_high']:.4f}), p={row['p']:.4g}"

    lines = [
        "# 全国节假日敏感性分析",
        "",
        "## 日历定义",
        "",
        f"使用 python-holidays {provenance['python_holidays_version']} 的 PUBLIC 类别与 `observed=True` 生成常规公共日期。巴西核心日历有 {provenance['date_counts']['python_holidays_core']['BR']} 日，墨西哥有 {provenance['date_counts']['python_holidays_core']['MX']} 日。再依据TSE及INE官方日历补入两国全国联邦选举日，最终巴西 {provenance['date_counts']['after_national_election_supplement']['BR']} 日、墨西哥 {provenance['date_counts']['after_national_election_supplement']['MX']} 日。该指标按公共日期或全国联邦选举日定义，用于医疗排程敏感性。",
        "",
        "该指标采用标准PUBLIC日期与两国全国联邦选举日的并集。巴西Good Friday日期由标准日历和年度联邦日历核对，其宗教假日法律框架涉及地方规定。州/市的其他节日、当地选举休假、行政机关自由裁量假日、调休桥接日及医院自身排班均不在指标范围内。日期、节名、库版本和政府规则来源见 `holiday_calendar_dates.csv` 与 `holiday_calendar_provenance.json`。",
        "",
        "## 本批拟合前固定的六个追加敏感性模型",
        "",
        "| 模型 | 主要估计 | 入模日数 / 拟合保留 N |",
        "|---|---|---:|",
    ]
    for country, label in (("brazil", "巴西"), ("mexico", "墨西哥 ERA5-Land")):
        state_key = f"{country}_strict_n20_holiday_lag07"
        state_sum = summary_lookup[state_key]
        lines.append(
            f"| {label}严格 N20，温度 lag 0–7 + 假日 | 温度累积 RR {fmt(state_key, 'cumulative_temperature_lag_0_7_per_C')}; 假日 RR {fmt(state_key, 'national_holiday_rate_ratio')} | {state_sum['input_rows']:,} / {state_sum['model_rows']:,} |"
        )
    for outcome in ("any_drainage", "any_removal"):
        model = f"sp_{outcome}_holiday_lag07"
        summary = summary_lookup[model]
        lines.append(
            f"| 巴西 SP {outcome}，温度与火 PM2.5 各 lag 0–7 + 假日 | 温度累积 RR {fmt(model, 'cumulative_temperature_lag_0_7_per_C')}; 火 PM2.5 累积 RR {fmt(model, 'cumulative_fire_PM25_lag_0_7_per_10_units')}; 假日 RR {fmt(model, 'national_holiday_rate_ratio')} | {summary['input_rows']:,} / {summary['model_rows']:,} |"
        )
    for phenotype, model in [
        ("广义 N20–N23", "broad_urgent_elective_holiday_joint"),
        ("严格 N20", "strict_urgent_elective_holiday_joint"),
    ]:
        summary = summary_lookup[model]
        lines.append(
            f"| {phenotype} 急诊/择期联合模型 | heat-only 急诊 RR / 择期 RR = {fmt(model, 'heat_only_urgent_effect_over_elective_effect')}; 烟-only = {fmt(model, 'smoke_only_urgent_effect_over_elective_effect')}; both = {fmt(model, 'both_urgent_effect_over_elective_effect')} | {summary['input_rows']:,} / {summary['model_rows']:,} |"
        )
    lines.extend([
        "",
        "急诊/择期联合模型同时估计 elective 与 urgent 的全国假日 RR；两端点假日效应及完整协方差线性对比见 `holiday_model_contrasts.csv` 和 `holiday_full_covariance.json`。状态模型按状态与日期双向 CRV1 聚类，累积温度与假日区间使用 t(min(G_state,G_date)−1)；其余模型按既有定义使用健康区域 CRV1。所有滞后项、假日项、分层样本、全部正负系数和完整协方差均已保存。",
        "",
        "这些模型检验全国性假日控制加入后既定结局规格的稳健性；它们不表示州/市节日、桥接日或医院实际运营安排已被完整控制。",
        "",
    ])
    (OUT / "holiday_report_zh.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    calendar_frame, provenance = build_holiday_calendar()
    specifications = {
        "study_window": "2018-01-01 through 2022-12-31",
        "six_fits_only": SPECIFICATIONS,
        "holiday_covariate": "recurring PUBLIC date or nationwide federal election; country-specific dates in holiday_calendar_dates.csv",
        "calendar_file": str(OUT / "holiday_calendar_dates.csv"),
        "calendar_provenance_file": str(OUT / "holiday_calendar_provenance.json"),
    }
    (OUT / "holiday_model_specifications.json").write_text(
        json.dumps(specifications, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8"
    )
    print("HOLIDAY_CALENDAR_LOCKED", json.dumps(provenance["date_counts"], ensure_ascii=True), flush=True)

    mexico_module = import_module("mexico_strict_phenotype_helper", WORK / "mexico_strict_phenotype.py")
    sp_module = import_module("sp_sensitivity_helper", WORK / "fit_sp_sensitivity.py")
    holiday_calendar_global[0] = calendar_frame

    all_coefficients, all_contrasts, all_summaries, all_models = [], [], [], {}

    def persist_fit_checkpoint(status="in_progress"):
        checkpoint = {
            "status": status,
            "models_executed": list(all_models),
            "coefficients": all_coefficients,
            "contrasts": all_contrasts,
            "summaries": all_summaries,
            "models": all_models,
        }
        (OUT / "holiday_fit_checkpoint.json").write_text(
            json.dumps(checkpoint, indent=2, ensure_ascii=False, allow_nan=False),
            encoding="utf-8",
        )

    state_panel = state_level_frame(mexico_module)
    state_contrasts = []
    for country in ("BR", "MX"):
        country_panel = state_panel.loc[state_panel.country.eq(country)].copy()
        state_coef, state_con, state_summary, state_details = fit_state_strict(
            country_panel, country, calendar_frame
        )
        all_coefficients.extend(state_coef)
        all_contrasts.extend(state_con)
        state_contrasts.extend(state_con)
        all_summaries.append(state_summary)
        all_models[state_summary["model"]] = state_details
        persist_fit_checkpoint()

    sp_full, sp_lagged, sp_lags, sp_audit = load_sp_frame(sp_module, calendar_frame)

    def collect_sp_fit(coef_rows, contrast_rows, summary, details):
        all_coefficients.extend(coef_rows)
        all_contrasts.extend(contrast_rows)
        all_summaries.append(summary)
        all_models[summary["model"]] = details
        persist_fit_checkpoint()

    _, sp_con, _, _ = fit_sp_sensitivities(
        sp_full, sp_lagged, sp_lags, sp_audit, on_fit=collect_sp_fit
    )

    brazil_joint = brazil_urgent_elective_panel(sp_full)
    broad_coef, broad_con, broad_summary, broad_details = fit_joint_urgent_elective(brazil_joint, "broad")
    all_coefficients.extend(broad_coef)
    all_contrasts.extend(broad_con)
    all_summaries.append(broad_summary)
    all_models[broad_summary["model"]] = broad_details
    persist_fit_checkpoint()

    strict_coef, strict_con, strict_summary, strict_details = fit_joint_urgent_elective(brazil_joint, "strict")
    all_coefficients.extend(strict_coef)
    all_contrasts.extend(strict_con)
    all_summaries.append(strict_summary)
    all_models[strict_summary["model"]] = strict_details
    persist_fit_checkpoint()

    if len(all_models) != 6:
        raise ValueError(f"Expected exactly six holiday fits, found {len(all_models)}")
    expected_models = {
        "brazil_strict_n20_holiday_lag07",
        "mexico_strict_n20_holiday_lag07",
        "sp_any_drainage_holiday_lag07",
        "sp_any_removal_holiday_lag07",
        "broad_urgent_elective_holiday_joint",
        "strict_urgent_elective_holiday_joint",
    }
    if set(all_models) != expected_models:
        raise ValueError(f"Executed model set differs from locked six fits: {sorted(all_models)}")

    pd.DataFrame(all_coefficients).to_csv(
        OUT / "holiday_model_coefficients.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(all_contrasts).to_csv(
        OUT / "holiday_model_contrasts.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(all_summaries).to_csv(
        OUT / "holiday_model_summary.csv", index=False, encoding="utf-8-sig"
    )

    sample_audit = []
    for country, state_days in state_panel.groupby("country", observed=True):
        relevant_holiday_dates = pd.to_datetime(
            calendar_frame.loc[calendar_frame.country.eq(country), "date"]
        )
        sample_audit.append({
            "analysis": f"{COUNTRY_NAMES[country]}_strict_n20_holiday_lag07",
            "country": country,
            "geographies": int(state_days.country_state.nunique()),
            "calendar_rows": int(len(state_days)),
            "unique_dates": int(state_days.date.nunique()),
            "events": int(state_days.strict_N20.sum()),
            "national_holiday_dates": int(relevant_holiday_dates.nunique()),
            "holiday_geography_days": int(state_days.date.isin(relevant_holiday_dates).sum()),
            "input_rows_after_lag_warmup": int(state_days.country_state.nunique() * (1826 - 7)),
        })
    for summary in all_summaries:
        sample_audit.append({
            "analysis": summary["model"],
            "country": summary.get("country", "BR"),
            "geographies": summary["cluster_counts"].get("hrid", summary["cluster_counts"].get("country_state")),
            "calendar_rows": summary["input_rows"],
            "unique_dates": 1826,
            "events": summary["input_outcome_events"],
            "national_holiday_dates": None,
            "holiday_geography_days": summary["holiday_input_rows"],
            "input_rows_after_lag_warmup": summary["input_rows"],
        })
    pd.DataFrame(sample_audit).to_csv(
        OUT / "holiday_sample_audit.csv", index=False, encoding="utf-8-sig"
    )

    manifest = {
        "script": str(Path(__file__).resolve()),
        "window": specifications["study_window"],
        "holiday_calendar": provenance,
        "model_specifications_locked_before_fit": SPECIFICATIONS,
        "models_executed": list(all_models),
        "state_input_sources": {
            "Brazil": str(PANEL_PATH),
            "Mexico_ERA5": str(mexico_module.MEXICO_ERA5_PATH),
            "clinical_aggregate": str(CLINICAL_PATH),
        },
        "SP_link_source": str(SP_LINK_PATH),
        "SP_link_sha256": sp_audit["SP_input_sha256"],
        "outputs": [
            "holiday_calendar_dates.csv",
            "holiday_calendar_provenance.json",
            "holiday_model_specifications.json",
            "holiday_model_coefficients.csv",
            "holiday_model_contrasts.csv",
            "holiday_model_summary.csv",
            "holiday_sample_audit.csv",
            "holiday_full_covariance.json",
            "holiday_run_manifest.json",
            "holiday_report_zh.md",
            "holiday_fit_checkpoint.json",
        ],
        "models": all_models,
    }
    (OUT / "holiday_full_covariance.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8"
    )
    (OUT / "holiday_run_manifest.json").write_text(
        json.dumps({key: value for key, value in manifest.items() if key != "models"},
                   indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8"
    )
    make_report(calendar_frame, state_contrasts, sp_con, broad_con + strict_con, all_summaries, provenance)
    persist_fit_checkpoint(status="complete")
    print("HOLIDAY_SENSITIVITY_COMPLETE", flush=True)
    print(pd.DataFrame(all_contrasts).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
