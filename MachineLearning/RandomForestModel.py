
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split,KFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, roc_auc_score


# 读取数据集
data = pd.read_csv('/data/AllDataToBuildModel.csv',index_col=0)
# data.drop(columns=['Row.names'], inplace=True)

# 将基因表达量作为特征矩阵X，将样本的分组情况作为目标变量y
X = data.iloc[:, :-1].values
y = data.iloc[:, -1].values

# 设置交叉验证参数
n_splits = 5
kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

# 定义随机森林模型参数
n_estimators = 100
max_depth = 10
min_samples_split = 2
min_samples_leaf = 1

# 进行交叉验证
auc_scores = []
for train_idx, test_idx in kf.split(X):
    # 将数据集分为训练集和测试集
    X_train, y_train = X[train_idx], y[train_idx]
    X_test, y_test = X[test_idx], y[test_idx]

    # 训练多分类随机森林模型
    rf = RandomForestClassifier(n_estimators=n_estimators, max_depth=max_depth, 
                                 min_samples_split=min_samples_split, min_samples_leaf=min_samples_leaf,
                                 random_state=42)
    rf.fit(X_train, y_train)

    # 使用测试集计算模型AUC值
    y_pred_prob = rf.predict_proba(X_test)
    auc = roc_auc_score(y_test, y_pred_prob, multi_class='ovo')
    auc_scores.append(auc)

# 输出交叉验证结果
print("交叉验证AUC值:", auc_scores)
print("平均AUC值:", np.mean(auc_scores))

# # 构建随机森林分类模型
classifier = RandomForestClassifier(n_estimators=n_estimators, max_depth=max_depth, 
                                 min_samples_split=min_samples_split, min_samples_leaf=min_samples_leaf,
                                 random_state=42)
classifier.fit(X, y)

# 计算每个变量基因的贡献度
importances = classifier.feature_importances_
importances_df = pd.DataFrame(importances, index=data.columns[:-1], columns=['importance'])
importances_df.sort_values(by='importance', ascending=False, inplace=True)
print("变量基因的贡献度：\n", importances_df.iloc[0:10,:])


