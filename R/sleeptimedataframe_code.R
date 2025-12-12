# #create test dataframe
# #time starts at hr 12 and goes to hr 36 (noon to noon)
# #Vector of 0 and 1 with 1 being sleep
#
# s_df <- data.frame("time"=seq(from = 12, to = 36, by = .1),
#                    "S" = 0)
#
# #convert to sleep starting at 10pm and switching to wake at 6:30am
# s_df$S[s_df$time >= 22 & s_df$time < 30.5] <- 1
#
# #Objective 1: return a list with 2 entries
# # 1) the first one indicates the duration of sleep in hours (numeric values only)
# # 2) a named vector with sleep onset, sleep midpoint, and sleep offset times
#
# s_df[order(s_df$time), ] #sort ascending
#
# sleep_first <- match("1", s_df$S) #find first instance of sleep
# sleep_on <- s_df[sleep_first, "time"] #time of sleep onset
#
# sleep_last <- max(which(s_df$S == 1)) #find last instance of sleep
# sleep_off <- s_df[sleep_last, "time"] + .1 #awake onset
#
# sleep_duration <- sleep_off - sleep_on #duration of sleep in hours
#
# sleep_mid <- sleep_on + (sleep_duration/2) #midpoint of sleep
#
# sleep_times <- c(Sleep_on = sleep_on, sleep_mid = sleep_mid, sleep_off = sleep_off)
# print(sleep_times) #return a named vector with sleep times
#
# sleeptimes_list <- list(
#   sleep_duration,
#   sleep_times
# )
# print(sleeptimes_list)
#
# ## Answers
# # Sleep duration = 8.5
# # Sleep times: these can be formatted in different ways, but something like:
# # c(sleep_on = 20, sleep_mid = 26.25, sleep_off = 30.5)
# ###KHB note: sleep_on=22 (10pm)?
#
# ###############################################################################
# ## Objective 2: different formats for sleep timing
# # See if you can convert the sleep times from objective one to the following formats:
# # c(sleep_on = 22, sleep_mid = 2.25, sleep_off = 6.5) # hint, use modulus division (%%)
#
# sleeptime_num <- c(sleep_on, sleep_mid %% 12, sleep_off %% 12) #relies on start being PM and mid/off being AM
# print(sleeptime_num)
#
# # c(sleep_on = 10:00 pm, sleep_mid = 2:15 am, sleep_off = 6:30 am) # hint, will need to convert to strings
#
# hours <- floor(sleeptime_num)
# minutes <- round((sleeptime_num - hours) * 60)
# string_24hr <- paste0(sprintf("%02d", hours), ":", sprintf("%02d", minutes))
#
# time_posixlt <- strptime(string_24hr, format = "%H:%M")
# string_12hr <- format(time_posixlt, format = "%I:%M %p")
# print(string_12hr)
#
#
# ################################################################################
# # Objective 3: Return sleep midpoint in radians instead of clock time
# # time-of-day can be thought of as a circular variable. and represented in
# # either degrees or radian. For example, the first midnight can be set to 0 degrees
# # (or 0 radians) and the second midnight is 360 degrees (or 2*pi radians).
# # See if you can convert sleep midpoint clock time (from Objective #1) to radians.
#
# sleeptime_mid <- sleep_mid %% 12
# period_hours <- 24
# proportion <- sleeptime_mid / period_hours
# radians_mid <-proportion * (2 * pi)
# print(radians_mid)
#
# #Answer = ~0.5890486
#
#
#
#
#
#
#
#
#
#
#
#
#





