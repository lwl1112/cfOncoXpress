import os
import json
import pickle
import optuna
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import spearmanr
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.svm import SVR
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from xgboost import XGBRegressor

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

# =========================
# 0. Settings
# =========================

OUTDIR = "optuna_leave_one_control_out"
os.makedirs(OUTDIR, exist_ok=True)

feature_cols = [
    'meancov1k',
    'meancovNDR',
    'meanFFT950',
    'meanFFTNDR',
    'meanSL1k',
    'meansimpson1k_100_300'
]

use_cols = ['gene_name'] + feature_cols + ['TPM']

N_TRIALS_CLASSICAL = 150
N_TRIALS_NN = 80
RANDOM_STATE = 42

# =========================
# 1. Load controls
# =========================

control01 = pd.read_csv(
    '/athena/khuranalab/scratch/wel4007/csv_files/finalfeatures/redo/control01bigmat_all_rmdupgenes.csv'
)[use_cols]

control02 = pd.read_csv(
    '/athena/khuranalab/scratch/wel4007/csv_files/finalfeatures/redo/control02bigmat_all_rmdupgenes.csv'
)[use_cols]

# IMPORTANT: check this filename
wcmfc01 = pd.read_csv(
    '/athena/khuranalab/scratch/wel4007/csv_files/finalfeatures/redo/WCMFC01bigmat_all.csv'
)[use_cols]

datasets = {
    "control01": control01,
    "control02": control02,
    "wcmfc01": wcmfc01
}

# =========================
# 2. Helper functions
# =========================

def prepare_xy(df):
    df = df.dropna(subset=feature_cols + ["TPM"]).copy()

    X = df[feature_cols].values
    y = np.log2(df["TPM"].values + 1)

    return X, y


def calc_metrics(y_true, y_pred):
    rho = spearmanr(y_true, y_pred).correlation
    if np.isnan(rho):
        rho = -1

    return {
        "spearman": rho,
        "mse": mean_squared_error(y_true, y_pred),
        "mae": mean_absolute_error(y_true, y_pred)
    }


def save_pareto_plot(study, model_name):
    pareto_trials = study.best_trials

    pareto_df = pd.DataFrame({
        "spearman": [t.values[0] for t in pareto_trials],
        "mse": [t.values[1] for t in pareto_trials]
    })

    pareto_df.to_csv(
        os.path.join(OUTDIR, f"{model_name}_pareto_trials.csv"),
        index=False
    )

    plt.figure(figsize=(6, 5))
    plt.scatter(pareto_df["mse"], pareto_df["spearman"])
    plt.xlabel("MSE")
    plt.ylabel("Spearman")
    plt.title(f"{model_name} Pareto Front")
    plt.tight_layout()
    plt.savefig(
        os.path.join(OUTDIR, f"{model_name}_pareto_front.png"),
        dpi=300
    )
    plt.close()


def save_best_params(study, model_name):
    best_trial = max(
        study.best_trials,
        key=lambda t: t.values[0]
    )

    with open(os.path.join(OUTDIR, f"{model_name}_best_params.json"), "w") as f:
        json.dump(best_trial.params, f, indent=4)

    return best_trial

# =========================
# 3. Classical model builders
# =========================

def build_gbr(trial):
    return GradientBoostingRegressor(
        n_estimators=trial.suggest_int("n_estimators", 100, 1000),
        learning_rate=trial.suggest_float("learning_rate", 0.01, 0.3, log=True),
        max_depth=trial.suggest_int("max_depth", 2, 6),
        min_samples_split=trial.suggest_int("min_samples_split", 2, 20),
        min_samples_leaf=trial.suggest_int("min_samples_leaf", 1, 20),
        subsample=trial.suggest_float("subsample", 0.5, 1.0),
        random_state=RANDOM_STATE
    )


def build_xgb(trial):
    return XGBRegressor(
        n_estimators=trial.suggest_int("n_estimators", 100, 1000),
        max_depth=trial.suggest_int("max_depth", 2, 8),
        learning_rate=trial.suggest_float("learning_rate", 0.01, 0.3, log=True),
        subsample=trial.suggest_float("subsample", 0.5, 1.0),
        colsample_bytree=trial.suggest_float("colsample_bytree", 0.5, 1.0),
        reg_alpha=trial.suggest_float("reg_alpha", 1e-8, 10.0, log=True),
        reg_lambda=trial.suggest_float("reg_lambda", 1e-8, 10.0, log=True),
        objective="reg:squarederror",
        random_state=RANDOM_STATE,
        n_jobs=-1
    )


def build_rf(trial):
    return RandomForestRegressor(
        n_estimators=trial.suggest_int("n_estimators", 100, 1000),
        max_depth=trial.suggest_int("max_depth", 2, 30),
        min_samples_split=trial.suggest_int("min_samples_split", 2, 20),
        min_samples_leaf=trial.suggest_int("min_samples_leaf", 1, 20),
        max_features=trial.suggest_categorical("max_features", ["sqrt", "log2", None]),
        random_state=RANDOM_STATE,
        n_jobs=-1
    )


def build_svr(trial):
    return Pipeline([
        ("scaler", StandardScaler()),
        ("model", SVR(
            kernel=trial.suggest_categorical("kernel", ["rbf", "linear", "poly"]),
            C=trial.suggest_float("C", 1e-2, 1e3, log=True),
            epsilon=trial.suggest_float("epsilon", 1e-3, 1.0, log=True),
            gamma=trial.suggest_categorical("gamma", ["scale", "auto"])
        ))
    ])


classical_builders = {
    "GBR": build_gbr,
    "XGB": build_xgb,
    "RF": build_rf,
    "SVR": build_svr
}

# =========================
# 4. Classical objective
# =========================

def make_classical_objective(build_fn):

    def objective(trial):
        spearman_list = []
        mse_list = []

        for heldout_name, test_df in datasets.items():

            train_df = pd.concat(
                [df for name, df in datasets.items() if name != heldout_name],
                ignore_index=True
            )

            X_train, y_train = prepare_xy(train_df)
            X_test, y_test = prepare_xy(test_df)

            model = build_fn(trial)
            model.fit(X_train, y_train)

            pred = model.predict(X_test)

            metrics = calc_metrics(y_test, pred)
            spearman_list.append(metrics["spearman"])
            mse_list.append(metrics["mse"])

        return np.mean(spearman_list), np.mean(mse_list)

    return objective

# =========================
# 5. Run classical studies
# =========================

all_best_trials = []

for model_name, build_fn in classical_builders.items():

    print(f"\nRunning Optuna for {model_name}")

    study = optuna.create_study(
        study_name=model_name,
        storage=f"sqlite:///{OUTDIR}/{model_name}.db",
        directions=["maximize", "minimize"],
        load_if_exists=True
    )

    study.optimize(
        make_classical_objective(build_fn),
        n_trials=N_TRIALS_CLASSICAL
    )

    trials_df = study.trials_dataframe()
    trials_df.to_csv(
        os.path.join(OUTDIR, f"{model_name}_trials.csv"),
        index=False
    )

    save_pareto_plot(study, model_name)
    best_trial = save_best_params(study, model_name)

    all_best_trials.append({
        "model": model_name,
        "best_spearman": best_trial.values[0],
        "best_mse": best_trial.values[1],
        "best_params": best_trial.params
    })

# =========================
# 6. NN model builders
# =========================

def build_basic_mlp(trial, input_dim):
    n_layers = trial.suggest_int("n_layers", 1, 3)
    lr = trial.suggest_float("learning_rate", 1e-4, 1e-2, log=True)

    model = keras.Sequential()
    model.add(layers.Input(shape=(input_dim,)))

    for i in range(n_layers):
        units = trial.suggest_categorical(f"units_l{i}", [16, 32, 64, 128])
        activation = trial.suggest_categorical(f"activation_l{i}", ["relu", "elu", "tanh"])

        model.add(layers.Dense(units, activation=activation))

        dropout = trial.suggest_float(f"dropout_l{i}", 0.0, 0.5)
        model.add(layers.Dropout(dropout))

    model.add(layers.Dense(1))

    model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=lr),
        loss="mse",
        metrics=["mae"]
    )

    return model


def build_regularized_mlp(trial, input_dim):
    n_layers = trial.suggest_int("n_layers", 1, 3)
    lr = trial.suggest_float("learning_rate", 1e-4, 1e-2, log=True)
    l2_strength = trial.suggest_float("l2", 1e-6, 1e-2, log=True)

    inputs = keras.Input(shape=(input_dim,))
    x = inputs

    for i in range(n_layers):
        units = trial.suggest_categorical(f"units_l{i}", [16, 32, 64, 128])
        activation = trial.suggest_categorical(f"activation_l{i}", ["relu", "elu", "tanh"])

        x = layers.Dense(
            units,
            activation=activation,
            kernel_regularizer=keras.regularizers.l2(l2_strength)
        )(x)

        use_bn = trial.suggest_categorical(f"batchnorm_l{i}", [True, False])
        if use_bn:
            x = layers.BatchNormalization()(x)

        dropout = trial.suggest_float(f"dropout_l{i}", 0.0, 0.5)
        x = layers.Dropout(dropout)(x)

    outputs = layers.Dense(1)(x)
    model = keras.Model(inputs, outputs)

    model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=lr),
        loss="mse",
        metrics=["mae"]
    )

    return model


def build_residual_mlp(trial, input_dim):
    width = trial.suggest_categorical("width", [32, 64, 128])
    n_blocks = trial.suggest_int("n_blocks", 1, 3)
    activation = trial.suggest_categorical("activation", ["relu", "elu", "tanh"])
    dropout = trial.suggest_float("dropout", 0.0, 0.5)
    l2_strength = trial.suggest_float("l2", 1e-6, 1e-2, log=True)
    lr = trial.suggest_float("learning_rate", 1e-4, 1e-2, log=True)

    inputs = keras.Input(shape=(input_dim,))

    x = layers.Dense(
        width,
        activation=activation,
        kernel_regularizer=keras.regularizers.l2(l2_strength)
    )(inputs)

    for b in range(n_blocks):
        shortcut = x

        x = layers.Dense(
            width,
            activation=activation,
            kernel_regularizer=keras.regularizers.l2(l2_strength)
        )(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(dropout)(x)

        x = layers.Dense(
            width,
            activation=None,
            kernel_regularizer=keras.regularizers.l2(l2_strength)
        )(x)
        x = layers.BatchNormalization()(x)

        x = layers.Add()([x, shortcut])
        x = layers.Activation(activation)(x)

    outputs = layers.Dense(1)(x)

    model = keras.Model(inputs, outputs)

    model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=lr),
        loss="mse",
        metrics=["mae"]
    )

    return model


nn_builders = {
    "Basic_MLP": build_basic_mlp,
    "Regularized_MLP": build_regularized_mlp,
    "Residual_MLP": build_residual_mlp
}

# =========================
# 7. NN objective
# =========================

def make_nn_objective(build_fn):

    def objective(trial):
        spearman_list = []
        mse_list = []

        batch_size = trial.suggest_categorical("batch_size", [16, 32, 64])

        for heldout_name, test_df in datasets.items():

            tf.keras.backend.clear_session()

            train_df = pd.concat(
                [df for name, df in datasets.items() if name != heldout_name],
                ignore_index=True
            )

            X_train, y_train = prepare_xy(train_df)
            X_test, y_test = prepare_xy(test_df)

            scaler = StandardScaler()
            X_train_scaled = scaler.fit_transform(X_train)
            X_test_scaled = scaler.transform(X_test)

            model = build_fn(trial, input_dim=X_train_scaled.shape[1])

            early_stop = keras.callbacks.EarlyStopping(
                monitor="val_loss",
                patience=10,
                restore_best_weights=True
            )

            model.fit(
                X_train_scaled,
                y_train,
                validation_split=0.2,
                epochs=200,
                batch_size=batch_size,
                verbose=0,
                callbacks=[early_stop]
            )

            pred = model.predict(X_test_scaled, verbose=0).flatten()

            metrics = calc_metrics(y_test, pred)
            spearman_list.append(metrics["spearman"])
            mse_list.append(metrics["mse"])

        return np.mean(spearman_list), np.mean(mse_list)

    return objective

# =========================
# 8. Run NN studies
# =========================

for model_name, build_fn in nn_builders.items():

    print(f"\nRunning Optuna for {model_name}")

    study = optuna.create_study(
        study_name=model_name,
        storage=f"sqlite:///{OUTDIR}/{model_name}.db",
        directions=["maximize", "minimize"],
        load_if_exists=True
    )

    study.optimize(
        make_nn_objective(build_fn),
        n_trials=N_TRIALS_NN
    )

    trials_df = study.trials_dataframe()
    trials_df.to_csv(
        os.path.join(OUTDIR, f"{model_name}_trials.csv"),
        index=False
    )

    save_pareto_plot(study, model_name)
    best_trial = save_best_params(study, model_name)

    all_best_trials.append({
        "model": model_name,
        "best_spearman": best_trial.values[0],
        "best_mse": best_trial.values[1],
        "best_params": best_trial.params
    })

# =========================
# 9. Compare all models
# =========================

best_summary_df = pd.DataFrame(all_best_trials)

best_summary_df = best_summary_df.sort_values(
    by=["best_spearman", "best_mse"],
    ascending=[False, True]
)

best_summary_df.to_csv(
    os.path.join(OUTDIR, "all_model_best_summary.csv"),
    index=False
)

# print(best_summary_df)

# =========================
# 10. Select final model
# =========================

best_row = best_summary_df.iloc[0]

final_model_name = best_row["model"]
final_best_params = best_row["best_params"]

print("Final selected model:", final_model_name)
print("Best Spearman:", best_row["best_spearman"])
print("Best MSE:", best_row["best_mse"])
print("Best params:", final_best_params)

# =========================
# 11. Refit final model on all controls
# =========================

all_controls = pd.concat(
    [control01, control02, wcmfc01],
    ignore_index=True
)

X_all, y_all = prepare_xy(all_controls)

fixed_trial = optuna.trial.FixedTrial(final_best_params)

if final_model_name in classical_builders:
    final_model = classical_builders[final_model_name](fixed_trial)
    final_model.fit(X_all, y_all)

    with open(os.path.join(OUTDIR, f"final_{final_model_name}.pkl"), "wb") as f:
        pickle.dump(final_model, f)

else:
    scaler = StandardScaler()
    X_all_scaled = scaler.fit_transform(X_all)

    final_model = nn_builders[final_model_name](
        fixed_trial,
        input_dim=X_all_scaled.shape[1]
    )

    batch_size = final_best_params.get("batch_size", 32)

    early_stop = keras.callbacks.EarlyStopping(
        monitor="loss",
        patience=10,
        restore_best_weights=True
    )

    final_model.fit(
        X_all_scaled,
        y_all,
        epochs=200,
        batch_size=batch_size,
        verbose=0,
        callbacks=[early_stop]
    )

    final_model.save(os.path.join(OUTDIR, f"final_{final_model_name}.keras"))

    with open(os.path.join(OUTDIR, f"final_{final_model_name}_scaler.pkl"), "wb") as f:
        pickle.dump(scaler, f)