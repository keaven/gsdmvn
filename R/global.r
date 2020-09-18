# These global variables are declared to eliminate associated R cmd check warnings.
# There is no other identified functional impact of these global declarations.
# These are column names used in input and/or output tibbles in the package

utils::globalVariables(
  c(
    'Events',
    'Probability',
    'Stratum',
    'Time',
    'AHR',
    'N',
    'Z',
    'dropoutRate',
    'failRate',
    'duration',
    'rate',
    'theta',
    'info',
    'info0',
    'Bound',
    'Analysis',
    'h'
  )
)
