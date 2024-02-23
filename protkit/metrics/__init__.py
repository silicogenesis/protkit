"""
Package `metrics` contains classes to perform various evaluations on biological data.

- `StructureEval` class in `structure_eval`: Various metrics such as RMSD, TMScore, GDT TS and Fnat scores to evaluate protein structures.
- `SequenceEval` class in `sequence_eval`: Various metrics to such as sequence identity, similarity, and coverage to evaluate protein sequences.
- `ClassificationEval` class in `classification_eval`: Various metrics such as accuracy, precision, recall, and F1 score to the results of protein classification problems.
- `RegressionEval` class in `regression_eval`: Various metrics such as mean squared error, mean absolute error, and R2 score to evaluate the results of protein regression problems.
- `UtilityEval` class in `utility_eval`: Various metrics such as memory usage, execution time, and disk usage to evaluate metrics related to the execution of a program.
"""

from .structure_eval import StructureEval
from .sequence_eval import SequenceEval
from .regression_eval import RegressionEval
from .classification_eval import ClassificationEval
from .utility_eval import UtilityEval

from .scoring_matrix import ScoringMatrix