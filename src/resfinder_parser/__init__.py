__version__ = "0.1.0"

from .data_classes import Phenotype, IsolatePhenotypes, SeqRegion, IsolateSummary
from .resfinder_result_parser import ResfinderParser, ResfinderCollector

__all__ = [
    "__version__",
    "Phenotype",
    "IsolatePhenotypes", 
    "SeqRegion",
    "IsolateSummary",
    "ResfinderParser",
    "ResfinderCollector",
]
