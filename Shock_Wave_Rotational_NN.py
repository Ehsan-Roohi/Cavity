# Import necessary libraries
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import matplotlib.pyplot as plt
import os
from sklearn.preprocessing import MinMaxScaler

# -- 1. Data Loading and Normalization --

# --- Path to the directory containing the data files ---
DATASET_PATH = './'

# Define training and testing files for interpolation
train_files = ['M14.txt', 'M15.txt', 'M16.txt', 'M18.txt', 'M19.txt', 'M20.txt']
test_file = 'M17.txt'
# ---------------------------------------------------------

# Check if data files exist
if not os.path.exists(os.path.join(DATASET_PATH, train_files[0])):
    raise FileNotFoundError(f"Data files not found in '{DATASET_PATH}'. Please ensure the .txt files are uploaded.")


def load_data(files, path):
    """
    This function loads data from text files.
    Column Indexing:
    - 0: X/MFP (Position)
    - 2: VELOCITY
    - 4: TRANS TEMP
    - 5: ROT TEMP
    """
    x_data, y_data, mach_data = [], [], []
    for file_name in files:
        mach_number = float(file_name.replace('M', '').replace('.txt', ''))
        full_path = os.path.join(path, file_name)
        data = np.loadtxt(full_path)
        x_data.extend(data[:, 0])
        mach_data.extend([mach_number] * len(data))
        y_data.extend(data[:, [2, 4, 5]])
        
    return np.array(x_data), np.array(y_data), np.array(mach_data)

print("Loading data...")
x_train, y_train_original, mach_train = load_data(train_files, DATASET_PATH)
x_test, y_test_original, mach_test = load_data([test_file], DATASET_PATH)

X_train_original = np.stack([x_train, mach_train], axis=1)
X_test_original = np.stack([x_test, mach_test], axis=1)

print("Normalizing data...")
x_scaler = MinMaxScaler()
y_scaler = MinMaxScaler()

X_train_scaled = x_scaler.fit_transform(X_train_original)
y_train_scaled = y_scaler.fit_transform(y_train_original)

X_test_scaled = x_scaler.transform(X_test_original)
print("Normalization successful.")


# -- 2. Neural Network Design with Dropout --
print("\nBuilding the optimized model for memory...")
model = keras.Sequential([
    layers.Dense(128, activation='relu', input_shape=[2]),
    layers.Dropout(0.05), # Dropout rate reduced to make the confidence interval narrower
    layers.Dense(128, activation='relu'),
    layers.Dropout(0.05), # Dropout rate reduced to make the confidence interval narrower
    layers.Dense(128, activation='relu'),
    layers.Dropout(0.05), # Dropout rate reduced to make the confidence interval narrower
    layers.Dense(3)
])

initial_learning_rate = 0.001
lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
    initial_learning_rate,
    decay_steps=20000,
    decay_rate=0.9,
    staircase=True)

model.compile(optimizer=keras.optimizers.Adam(learning_rate=lr_schedule), 
              loss='mean_squared_error')
model.summary()

print("\nStarting training for 3000 epochs with learning rate decay...")
history = model.fit(
    X_train_scaled, y_train_scaled,
    epochs=3000,
    validation_split=0.2,
    verbose=1,
    batch_size=16
)
print("Training finished.")


# -- 3. Evaluation, Uncertainty Estimation, and Plotting --
print("\nEvaluating and estimating uncertainty with MC Dropout...")

# Increased number of passes for a smoother, more stable interval estimate
n_passes = 200
predictions_mc = []
for i in range(n_passes):
    print(f"  MC Pass {i+1}/{n_passes}")
    # By setting training=True, we keep dropout active during prediction
    predictions_mc.append(model(X_test_scaled, training=True))

# Stack predictions to shape (n_passes, n_samples, n_outputs)
predictions_mc = tf.stack(predictions_mc)

# Calculate mean and standard deviation across the passes
mean_prediction_scaled = tf.reduce_mean(predictions_mc, axis=0)
std_dev_scaled = tf.math.reduce_std(predictions_mc, axis=0)

# Inverse transform the mean prediction to the original scale
y_pred_mean = y_scaler.inverse_transform(mean_prediction_scaled)

# Calculate 95% confidence interval bounds in the scaled domain and then transform
lower_bound_scaled = mean_prediction_scaled - 1.96 * std_dev_scaled
upper_bound_scaled = mean_prediction_scaled + 1.96 * std_dev_scaled
y_pred_lower = y_scaler.inverse_transform(lower_bound_scaled)
y_pred_upper = y_scaler.inverse_transform(upper_bound_scaled)

print("Saving plots...")
test_mach_number = mach_test[0]
font_settings = {'fontsize': 16, 'fontweight': 'bold'}
label_font_settings = {'fontsize': 14}

# Sort the test data for correct plotting of confidence interval
sort_indices = np.argsort(x_test)
x_test_sorted = x_test[sort_indices]
y_test_sorted = y_test_original[sort_indices]
y_pred_mean_sorted = y_pred_mean[sort_indices]
y_pred_lower_sorted = y_pred_lower[sort_indices]
y_pred_upper_sorted = y_pred_upper[sort_indices]

# --- Loss Curve Plot (Logarithmic Scale) ---
plt.figure(figsize=(12, 8))
plt.plot(history.history['loss'], label='Training Loss')
plt.plot(history.history['val_loss'], label='Validation Loss')
plt.yscale('log')
plt.title('Model Training & Validation Loss (Log Scale)', **font_settings)
plt.xlabel('Epoch', **label_font_settings)
plt.ylabel('Loss (Mean Squared Error)', **label_font_settings)
plt.legend(fontsize=12)
plt.grid(True, which="both", ls="--")
plt.tight_layout()
plt.savefig('loss_curve_log.jpg', dpi=300)
plt.close()

# --- Comparison Plot for Velocity ---
plt.figure(figsize=(12, 8))
plt.scatter(x_test_sorted, y_test_sorted[:, 0], label='DSMC Solution', color='blue', s=20)
plt.plot(x_test_sorted, y_pred_mean_sorted[:, 0], label='Mean Prediction', color='red', linewidth=2)
plt.fill_between(x_test_sorted, y_pred_lower_sorted[:, 0], y_pred_upper_sorted[:, 0], color='red', alpha=0.2, label='95% Confidence Interval')
plt.title(f'Velocity Profile Comparison (Mach = {test_mach_number})', **font_settings)
plt.xlabel('Position (X/MFP)', **label_font_settings)
plt.ylabel('Normalized Velocity', **label_font_settings)
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
plt.savefig('comparison_velocity.jpg', dpi=300)
plt.close()

# --- Comparison Plot for Translational Temperature ---
plt.figure(figsize=(12, 8))
plt.scatter(x_test_sorted, y_test_sorted[:, 1], label='DSMC Solution', color='blue', s=20)
plt.plot(x_test_sorted, y_pred_mean_sorted[:, 1], label='Mean Prediction', color='green', linewidth=2)
plt.fill_between(x_test_sorted, y_pred_lower_sorted[:, 1], y_pred_upper_sorted[:, 1], color='green', alpha=0.2, label='95% Confidence Interval')
plt.title(f'Translational Temperature Profile (Mach = {test_mach_number})', **font_settings)
plt.xlabel('Position (X/MFP)', **label_font_settings)
plt.ylabel('Translational Temperature', **label_font_settings)
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
plt.savefig('comparison_trans_temp.jpg', dpi=300)
plt.close()

# --- Comparison Plot for Rotational Temperature ---
plt.figure(figsize=(12, 8))
plt.scatter(x_test_sorted, y_test_sorted[:, 2], label='DSMC Solution', color='blue', s=20)
plt.plot(x_test_sorted, y_pred_mean_sorted[:, 2], label='Mean Prediction', color='purple', linewidth=2)
plt.fill_between(x_test_sorted, y_pred_lower_sorted[:, 2], y_pred_upper_sorted[:, 2], color='purple', alpha=0.2, label='95% Confidence Interval')
plt.title(f'Rotational Temperature Profile (Mach = {test_mach_number})', **font_settings)
plt.xlabel('Position (X/MFP)', **label_font_settings)
plt.ylabel('Rotational Temperature', **label_font_settings)
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
plt.savefig('comparison_rot_temp.jpg', dpi=300)
plt.close()

print("All plots saved successfully.")