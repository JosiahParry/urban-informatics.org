#' library(sf)
#' library(spdep)
#' library(tidyverse)
#' 
#' #' Spatial autocorrelation process:
#' #' 
#' #' - create neighbors
#' #' - create weights
#' #' - create lag
#' #' - Create moran plot
#' #'   - this is original variabl vs spatial lag
#' #' - calculate morans I for global
#' #' - local morans I for local autocorrelation
#' 
#' # Create  some data to work with
#' acs <- select(uitk::acs_raw, 
#'               fips = ct_id_10, med_house_income, 
#'               less_than_hs, by_pub_trans, bach) %>% 
#'   mutate(fips = as.character(fips),
#'          #median inpute for missing values
#'          med_house_income = replace_na(med_house_income, median(med_house_income, na.rm = TRUE)),
#'          less_than_hs = replace_na(less_than_hs, median(less_than_hs, na.rm = TRUE)),
#'          bach = replace_na(bach, median(bach, na.rm = TRUE)),
#'          by_pub_trans = replace_na(by_pub_trans, median(by_pub_trans, na.rm = TRUE))) 
#'   #filter(!is.na(med_house_income))
#' 
#' acs_sf <- left_join(uitk::suffolk_county, acs)
#' 
#' # is there a spatial relationship in median household income
#' # in boston?
#' 
#' ggplot(acs_sf, aes(fill = med_house_income)) +
#'   geom_sf(color = "black", lwd = 0.25) +
#'   theme_void() +
#'   labs(title = "Median Household Income in Boston")
#' 
#' 
#' 
#' # neighbors, weights, & lag -----------------------------------------------
#' 
#' acs_nb <- poly2nb(acs_sf)
#' 
#' acs_w <- nb2listw(acs_nb)
#' 
#' acs_lags <- acs_sf %>% 
#'   # create spatial lags
#'   mutate(across(where(is.numeric), ~lag.listw(acs_w, .x), .names = "{.col}_lag")) 
#' 
#' 
#' # moran plot --------------------------------------------------------------
#' #' Looks at the relationship between the spatial lag and the observed variable
#' #' The moran plot checks to see if observations are locally or globally 
#' #' different from the observed??? (check this)
#' #' 
#' 
#' # moran plot
#' acs_lags %>% 
#'   mutate(y = scale(med_house_income), 
#'          x = scale(med_house_income_lag),
#'          quadrant = case_when(
#'            x > 0 & y > 0 ~ "I",
#'            x > 0 & y < 0 ~ "II",
#'            x < 0 & y < 0 ~ "III",
#'            x < 0 & y > 0 ~ "IV"
#'          )) %>% 
#'   ggplot(aes(x, y, color = quadrant)) +
#'   geom_point() +
#'   theme_light() +
#'   geom_vline(xintercept = 0) +
#'   geom_hline(yintercept = 0) +
#'   labs(y = "Spatially lagged median income",
#'        x = "Observed median income",
#'        title = "Median HH income Moran Plot",
#'        caption = "One point removed to improve visualization extents.") +
#'   scale_color_manual(values = cpcinema::available_pals$joker[c(2,3,5,7)]) +
#'   lims(y = c(-1,3)) 
#' 
#' 
#' 
#' 
#' # density plot
#' acs_lags %>% 
#'   mutate(med_inc_lag_std = scale(med_house_income_lag)) %>% 
#'   ggplot(aes(med_inc_lag_std)) +
#'   geom_density() +
#'   geom_vline(xintercept = 0, color = "red", lty = 3) +
#'   theme_light() 
#' 
#' 
#' # Moran's I - global autocorrelation --------------------------------------
#' #moran(acs_lags$med_house_income_lag, acs_w, n = nrow(acs_lags), S0 = nrow(acs_lags))
#' 
#' moran.test(acs_lags$med_house_income_lag, acs_w)
#' #' We have positive spatial auto correlation. we can take this to mean that things that are close together are more closely related than those that are distance. 
#' #' In this case it means that high med hh income is located close to high median house hold income
#' 
#' 
#' # LISA - local indicators of spatial association --------------------------
#' 
#' #localmoran(acs_lags$med_house_income_lag, acs_w)
#' 
#' acs_lisa <- acs_lags %>% 
#'   mutate(lisas = as.data.frame(localmoran(med_house_income_lag, acs_w))) %>% 
#'   as_tibble() %>% 
#'   unpack(lisas) %>% 
#'   janitor::clean_names() %>% 
#'   st_as_sf()
#' 
#' ggplot(acs_lisa, aes(ii)) +
#'   geom_density() +
#'   labs(title = "LISAs density plot") +
#'   theme_light() +
#'   xlim(c(-1.5, 4.5))
#' 
#' 
#' acs_lisa %>% 
#'   ggplot(aes(fill = ii)) +
#'   geom_sf(color = "black", lwd = 0.2) +
#'   scale_fill_binned(n.breaks = 5) +
#'   theme_minimal() 
#' 
#' 
#' 
#' acs_lisa %>% 
#'   mutate(i_q = ntile(ii, 10)) %>% 
#'   ggplot(aes(fill = ii)) +
#'   geom_sf(color = "black", lwd = 0.2) +
#'   scale_fill_binned(n.breaks = 10, low = "#0C1F49", high = "#ECFFFB") +
#'   theme_minimal() 
#' 
#' 
#' 
#' readable_p <- function(p) {
#'   case_when(
#'     p < 0.001 ~ "p < 0.001",
#'     p < 0.01 ~ "p < 0.01", 
#'     p <= 0.05 ~ "p <= 0.05",
#'     p <= 0.1 ~ "p <= 0.10",
#'     TRUE ~ "non-significant"
#'   )
#' }
#' 
#' 
#' #' lisa is the moran I for the neighbor list
#' 
#' 
#' readable_p(acs_lisa$pr_z_0)
#' 
#' # steps to calculate the lag for the first row
#' acs_w$weights[[1]]
#' acs_w$neighbours[[1]]
#' acs_sf$med_house_income[acs_w$neighbours[[1]]]
#' acs_w$weights[[1]] * acs_sf$med_house_income[acs_w$neighbours[[1]]]
#' sum(acs_w$weights[[1]] * acs_sf$med_house_income[acs_w$neighbours[[1]]])
#' 
#' 
#' # To look at spatial autocorrelation
#' # we need to look at each polygon and then look at its neighbors
#' # how do we define the neighbors? With polygons there are two key ways that this is done
#' # we are looking at adjacency, also called contiguity
#' # contiguity : the state of bordering or being in direct contact with something.
#' # Contiguity is named after chess pieces. If you've seen the queens gtambit you should know how these pieces move
#' # there is rook and queen contiguity. 
#' # rook contiguity checks to see if polygons touch horizontally and vertically
#' # queen looks for horizonal, diagonal, and vertical. Anything that touches it. 
#' # this is in reality how neighborhoods are.
#' 
#' #' we want to capture what is called the spatial lag.
#' #' The spatial lag "he spatial lag captures the behavior of a variable in the immediate surroundings of each location" 
#' #' https://geographicdata.science/book/notebooks/06_spatial_autocorrelation.html
#' #' 
#' # we want to create a list of neighbors to each place and find the average values between them. This is called the lag
#' 
#' # create a list of neighborig polygons based on queen contiguity
#' acs_nb <- poly2nb(acs_sf)
#' 
#' acs_w <- nb2listw(acs_nb)
#' 
#' inc_lag <- acs_sf %>% 
#'   mutate(inc_lag = lag.listw(acs_w, med_house_income)) %>% 
#'   select(-geometry, everything())
#' 
#' # https://geographicdata.science/book/notebooks/07_local_autocorrelation.html
#' 
#' 
#' # the lagged variables is just the sum(y * neighbor inputs)
#' sum(acs_w$weights[[1]] * acs_sf$med_house_income[acs_w$neighbours[[1]]])
#' 
#' 
#' inc_lag %>% 
#'   mutate(inc_lag_std =  scale(inc_lag), 
#'          med_house_inc_std = scale(med_house_income)) %>% 
#'   ggplot(aes(med_house_inc_std, inc_lag_std)) +
#'   geom_point() +
#'   geom_hline(yintercept = 0) +
#'   geom_vline(xintercept = 0) +
#'   theme_light()
#' 
#' # lisa (local moran)
#' 
#' 
#' iris %>% 
#'   mutate(x = janitor::clean_names(iris)) %>% 
#'   str()
#' 
#' 
#' lisa <- janitor::clean_names(as.data.frame(localmoran(inc_lag$inc_lag, acs_w)))
#' 
#' lisa <- inc_lag %>% 
#'   mutate(lisa = as.data.frame(localmoran(inc_lag, acs_w))) %>% 
#'   as_tibble() %>% 
#'   tidyr::unpack(lisa) %>% 
#'   janitor::clean_names() %>% 
#'   st_as_sf()
#' 
#' lisa %>% 
#'  # mutate(i = ntile(ii, 5)) %>% 
#'   ggplot(aes(fill = ii)) +
#'   geom_sf(color = "black", lwd = 0.2) +
#'   scale_fill_binned(n.breaks = 5) +
#'   theme_minimal()
#' 
#' # Spatial Lag Regression --------------------------------------------------
#' # spatially lagged dependant variables let us take into consideration not just the observed variable but also the surroundings too. 
#' # "First, one could think of 𝛾(lag) as simply the effect of a unit change in your average surroundings. 
#' # Changes in observation i has spillover effects, which affects neighbors / spatially lagged vars
#' 
#' # Regular OLS
#' reg_lm <- lm(med_house_income ~ bach + by_pub_trans, data = inc_lag)
#' 
#' # Independent lagged OLS
#' lag_y_lm <- lm(inc_lag ~ bach + by_pub_trans, data = inc_lag)
#' 
#' # satially lagged x models
#' # Dependent lagged OLS 
#' dep_lagged_sf <- inc_lag %>% 
#'   mutate(bach_lag = lag.listw(acs_w, bach),
#'          trans_lag = lag.listw(acs_w, by_pub_trans))
#' 
#' lag_x_lm <- lm(med_house_income ~ bach_lag + trans_lag, data = dep_lagged_sf)
#' 
#' # All lagged OLS
#' lag_xy_lm <- lm(inc_lag ~ bach_lag + trans_lag, data = dep_lagged_sf)
#' 
#' summary(lag_xy_lm)
#' 
#' 
#' # Using lmSLX from spatialreg
#' 
#' 
#' spatreg_lm <- lmSLX(med_house_income ~ bach + by_pub_trans, data = dep_lagged_sf, listw = acs_w)
#' 
#' summary(spatreg_lm)  
#' 
#' # recreate spatialreg
#' 
#' ols_lm <- lm(log(med_house_income) ~ bach*100 + bach_lag*100 + by_pub_trans*100 + trans_lag*100, data = dep_lagged_sf)
#' all_preds_lagged_y <- lm(inc_lag ~ bach + bach_lag + by_pub_trans + trans_lag, data = dep_lagged_sf)
#' 
#' summary(ols_lm)
#' summary(all_preds_lagged_y)
#' 
#' 
#' 
#' # spatial lag model -------------------------------------------------------
#' sar_lm <- spautolm(inc_lag ~ bach + by_pub_trans, 
#'                    data = dep_lagged_sf, listw = acs_w)
#' 
#' summary(sar_lm)
#' 
#' 
#' # error model
#' error_lm <- errorsarlm(log(inc_lag) ~ bach + by_pub_trans, 
#'          data = dep_lagged_sf, listw = acs_w,
#'          # durbin fits lags
#'          Durbin = TRUE)
#' 
#' summary(error_lm)
#' xs
#' 
#' 
#' 
#'
# ## TO cover in other post:
# 
# - spatial lag
# - local autocorrelation
# - regression
# 
# ------