na_omit <- function(par, data)
{
  # Show warnings
  if(any(is.na(data$Aoto)))
    warning("omitting NA values in data$Aoto")
  if(any(is.na(data$Loto)))
    warning("omitting NA values in data$Loto")
  if(any(is.na(data$Lrel)))
    warning("omitting NA values in data$Lrel")
  if(any(is.na(data$Lrec)))
    warning("omitting NA values in data$Lrec")
  if(any(is.na(data$liberty)))
    warning("omitting NA values in data$liberty")

  # Omit NA values in otolith data
  include_otoliths <- !is.na(data$Aoto) & !is.na(data$Loto)
  data$Aoto <- data$Aoto[include_otoliths]
  data$Loto <- data$Loto[include_otoliths]

  # Omit NA values in tagging data
  include_tags <- !is.na(data$Lrel) & !is.na(data$Lrec) & !is.na(data$liberty)
  data$Lrel <- data$Lrel[include_tags]
  data$Lrec <- data$Lrec[include_tags]
  data$liberty <- data$liberty[include_tags]
  par$log_age <- par$log_age[include_tags]

  list(par=par, data=data)
}
