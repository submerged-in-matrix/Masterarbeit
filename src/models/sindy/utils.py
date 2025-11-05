from pathlib import Path
import numpy as np
import pandas as pd

def load_voltage_dataset(path: str | Path) -> pd.DataFrame:
    raise NotImplementedError

def preprocess_voltage_df(df: pd.DataFrame):
    raise NotImplementedError

def save_equation_latex(latex_str: str, out_path: str | Path) -> None:
    raise NotImplementedError
