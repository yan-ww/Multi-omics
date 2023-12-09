
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, GridSearchCV, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score, roc_curve

# 1. 导入数据集
data = pd.read_csv('/data/TraitsToBuildModel.csv',index_col=0)
# data.drop(columns=['Row.names'], inplace=True)

# 2. 数据预处理
X = data.iloc[:, :-1].values
y = data.iloc[:, -1].values
scaler = StandardScaler()
X = scaler.fit_transform(X)

n_splits = 5
kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
auc_scores = []
accuracy_scores = []
for train_idx, test_idx in kf.split(X):
    # 将数据集分为训练集和测试集
    X_train, y_train = X[train_idx], y[train_idx]
    X_test, y_test = X[test_idx], y[test_idx]

    # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # 3. 构建SVM模型
    param_grid = {'C': [0.1, 1, 10, 100], 'gamma': [0.1, 1, 10, 100], 'kernel': ['rbf']}
    svm = SVC(decision_function_shape='ovo')
    svm_cv = GridSearchCV(svm, param_grid, cv=5)
    svm_cv.fit(X_train, y_train)
    best_params = svm_cv.best_params_
    svm = SVC(C=best_params['C'], gamma=best_params['gamma'], kernel=best_params['kernel'], decision_function_shape='ovo')

    # 4. 训练模型并进行预测
    svm.fit(X_train, y_train)
    y_pred = svm.predict(X_test)

    # 5. 计算模型的准确性和模型中各变量基因的贡献度
    accuracy = accuracy_score(y_test, y_pred)
    print('Accuracy:', accuracy)
    auc = roc_auc_score(y_test, y_pred,multi_class='ovo')
    print('AUC:', auc)
    auc_scores.append(auc)
    accuracy_scores.append(accuracy)

print("交叉验证AUC值:", auc_scores)
print("平均AUC值:", np.mean(auc_scores))
print("交叉验证Accuracy值:", accuracy_scores)
print("平均Accuracy值:", np.mean(accuracy_scores))
# coef = svm.dual_coef_
# print('Gene importance rank:', np.argsort(coef)[::-1])
