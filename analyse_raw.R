## Symptom onset to peak PFU/RNA
df_sym <- readRDS("symptoms.RDS")

sym_ids <- df_sym$id_sub[df_sym$include_sym]

c_pfu <-
  data.frame("id_sub" = character(),
             "max_day" = integer())
c_copy <-
  data.frame("id_sub" = character(),
             "max_day" = integer())
for (i in 1:N){
  pfu <- extract_data_pfu(i)$vl
  copy <- extract_data(i)$vl
  id <- unique(df_merged$id_sub[df_merged$obs_id == i])
  if (any(pfu > 0, na.rm = TRUE)) {
    c_pfu <- rbind(c_pfu,
                   data.frame("id_sub" = id,
                              "max_day" = which(pfu == max(pfu, na.rm =
                                                             TRUE))[1]))
  }
  if (any(copy>0, na.rm=TRUE)){
    c_copy <- rbind(c_copy,
                    data.frame("id_sub" = id,
                               "max_day" = which(copy == max(copy, na.rm =
                                                               TRUE))[1]))
  }
}

pfu_merge <- merge(c_pfu, df_sym, by = "id_sub") %>%
  filter(include_sym) %>%
  filter(!is.na(as.double(contact_onset))) %>%
  mutate(sym_to_peak = max_day - as.double(contact_onset))

quantile(pfu_merge$sym_to_peak, prob = 0.5)
quantile(pfu_merge$sym_to_peak, prob = 0.25)
quantile(pfu_merge$sym_to_peak, prob = 0.75)

write.csv(pfu_merge$id_sub, "PFU_ids.csv", row.names = FALSE)
# PFU 3 (3,5) N=35

copy_merge <- merge(c_copy, df_sym, by = "id_sub") %>%
  filter(include_sym) %>%
  filter(!is.na(as.double(contact_onset))) %>%
  mutate(sym_to_peak = max_day - as.double(contact_onset))

write.csv(copy_merge$id_sub, "copy_ids.csv", row.names = FALSE)


quantile(copy_merge$sym_to_peak, prob = 0.5, type = 3)
quantile(copy_merge$sym_to_peak, prob = 0.25, type=3)
quantile(copy_merge$sym_to_peak, prob = 0.75, type=3)

# Copy 3 (3, 6)

df_vacc <- data.frame("id_sub"=unique(df_merged$id_sub))
vacc_status <- c()
for (id in df_vacc$id_sub){
  vacc_status <- c(vacc_status, unique(df_merged$vaccinated[df_merged$id_sub==id]))
}
df_vacc$vacc_status <- vacc_status

copy_merge <- merge(c_copy, df_sym, by = "id_sub")
copy_merge <- merge(copy_merge, df_vacc, by = "id_sub") %>%
  filter(include_sym) %>%
  filter(!vacc_status) %>%
  filter(!id_sub %in% c("ATA0364","ATA0386")) %>%
  mutate(sym_to_peak = max_day - as.double(contact_onset))

quantile(copy_merge$sym_to_peak, prob = 0.5, type = 3)
quantile(copy_merge$sym_to_peak, prob = 0.25, type=3)
quantile(copy_merge$sym_to_peak, prob = 0.75, type=3)

# Vacc 3 (3, 4) N = 17
# Unvacc 3 (2, 6) N = 15

pfu_merge <- merge(c_pfu, df_sym, by = "id_sub")
pfu_merge <- merge(pfu_merge, df_vacc, by = "id_sub") %>%
  filter(include_sym) %>%
  filter(!vacc_status) %>%
  filter(!id_sub %in% c("ATA0364","ATA0386")) %>%
  mutate(sym_to_peak = max_day - as.double(contact_onset))

quantile(pfu_merge$sym_to_peak, prob = 0.5, type = 3)
quantile(pfu_merge$sym_to_peak, prob = 0.25, type=3)
quantile(pfu_merge$sym_to_peak, prob = 0.75, type=3)

# Vacc 3 (3, 4) N = 17
# Unvacc 3 (2, 6) N = 13

## Peak PFU to peak RNA
# (Execute above)

df_pfu_copy <- merge(x=c_copy, y=c_pfu, by="id_sub") %>%
  mutate(same_day = max_day.x == max_day.y) %>%
  mutate(next_day = max_day.x == max_day.y-1) %>%
  mutate(prev_day = max_day.x == max_day.y+1) %>%
  mutate(day_diff = max_day.x - max_day.y)


sum(df_pfu_copy$same_day)
sum(df_pfu_copy$next_day)
sum(df_pfu_copy$prev_day)

write.csv(select(df_pfu_copy, id_sub, day_diff), "PFU_Copy_peak.csv", row.names = FALSE)

# length of infectious window raw data
remove_raw <-
  unique(df_merged$obs_id[df_merged$id_sub %in% c(
    "ATA0075",
    "BUZ091",
    "BUZ093",
    "BUZ099",
    "FUZR242",
    "FUZR155",
    "BUZ129",
    "FUZR114",
    "ATA0207"
  )]) 

exclude_pfu_raw <- c(12,18,23,25,41,56)  

group_1_pfu_raw <- group_1[!group_1%in%c(exclude_pfu_raw, remove_raw)]
group_2_pfu_raw <- group_2[!group_2%in%c(exclude_pfu_raw, remove_raw)]
group_3_pfu_raw <- group_3[!group_3%in%c(exclude_pfu_raw, remove_raw)]
group_4_pfu_raw <- group_4[!group_4%in%c(exclude_pfu_raw, remove_raw)]


# c(group_1_pfu_raw, group_2_pfu_raw, group_3_pfu_raw)
df_raw_unvacc <-
  data.frame("id_sub" = character(), "win_length" = integer())
for (i in c(group_1_pfu_raw)) {
  pfu_i <- df_pfu$pfu[df_pfu$obs_id == i]
  # print(pfu_i[!is.na(pfu_i)])
  # first day
  pfu_first <- which(pfu_i > 1)[1]
  pfu_last <- tail(which(pfu_i > 1), 1)
  df_raw_unvacc <-
    rbind(df_raw_unvacc,
          data.frame(
            "id_sub" = unique(df_pfu$id_sub[df_pfu$obs_id == i]),
            "win_length" = pfu_last - pfu_first + 1
          ))
}

df_raw_vacc <-
  data.frame("id_sub" = character(), "win_length" = integer())
for (i in group_2_pfu_raw) {
  pfu_i <- df_pfu$pfu[df_pfu$obs_id == i]
  # print(pfu_i[!is.na(pfu_i)])
  # first day
  pfu_first <- which(pfu_i > 1)[1]
  pfu_last <- tail(which(pfu_i > 1), 1)
  df_raw_vacc <-
    rbind(df_raw_vacc,
          data.frame(
            "id_sub" = unique(df_pfu$id_sub[df_pfu$obs_id == i]),
            "win_length" = pfu_last - pfu_first + 1
          ))
}

df_raw_all <-
  data.frame("id_sub" = character(), "win_length" = integer())
for (i in c(group_1_pfu_raw,
            group_2_pfu_raw,
            group_3_pfu_raw,
            group_4_pfu_raw)) {
  pfu_i <- df_pfu$pfu[df_pfu$obs_id == i]
  # print(pfu_i[!is.na(pfu_i)])
  # first day
  pfu_first <- which(pfu_i > 1)[1]
  pfu_last <- tail(which(pfu_i > 1), 1)
  df_raw_all <-
    rbind(df_raw_all,
          data.frame(
            "id_sub" = unique(df_pfu$id_sub[df_pfu$obs_id == i]),
            "win_length" = pfu_last - pfu_first + 1
          ))
}

write.csv(df_raw_vacc$id_sub, "pfu_vaccinated.csv", row.names = FALSE)
write.csv(df_raw_unvacc$id_sub, "pfu_unvaccinated.csv", row.names = FALSE)


df_raw_unvacc$win_length
df_raw_vacc$win_length

wilcox.test(df_raw_unvacc$win_length,
            df_raw_vacc$win_length)

quantile(df_raw_all$win_length, prob=0.5)
quantile(df_raw_all$win_length, prob=0.25)
quantile(df_raw_all$win_length, prob=0.75)

quantile(df_raw$win_length, prob=0.5)
quantile(df_raw$win_length, prob=0.25)
quantile(df_raw$win_length, prob=0.75)

# unvaccinated
# 5 (3, 7) N=24

# Vaccinated
# 5(3, 7) N=18

# All
# 5 (3,7) N = 41

# write.csv(df_raw$id_sub, "raw_pfu_vaccinated.csv", row.names = FALSE)

# Infectious more than five days post PFU peak

exclude_pfu_raw <- c(12,18,23,25,41,56)  

id_seq <- seq(1,57)
id_seq <- id_seq[!id_seq%in%exclude_pfu_raw]

# c(group_1_pfu_raw, group_2_pfu_raw, group_3_pfu_raw)
df_raw_post <-
  data.frame("id_sub" = character(), "win_length" = integer())
for (i in id_seq) {
  pfu_i <- df_pfu$pfu[df_pfu$obs_id == i]
  # print(pfu_i[!is.na(pfu_i)])
  # first day
  pfu_max <- which(pfu_i == max(pfu_i, na.rm = TRUE))[1]
  pfu_last <- tail(which(pfu_i > 1), 1)
  df_raw_post <-
    rbind(df_raw_post,
          data.frame(
            "id_sub" = unique(df_pfu$id_sub[df_pfu$obs_id == i]),
            "win_length" = pfu_last - pfu_max
          ))
}

df_raw_post[df_raw_post$win_length>5,]
df_raw_post$win_length[df_raw_post$win_length>5]

# length of infectious window raw data, from symptom onset
remove_raw_decline <-
  unique(df_merged$obs_id[df_merged$id_sub %in% c(
    "ATA0075",
    "BUZ091",
    "BUZ093",
    "BUZ099",
    "FUZR242",
    "FUZR155",
    "BUZ129",
    "FUZHH21"
  )]) 

exclude_pfu_raw <- c(12,18,23,25,41,56)  

group_1_pfu_raw <- group_1[!group_1%in%c(exclude_pfu_raw, remove_raw_decline)]
group_2_pfu_raw <- group_2[!group_2%in%c(exclude_pfu_raw, remove_raw_decline)]
group_3_pfu_raw <- group_3[!group_3%in%c(exclude_pfu_raw, remove_raw_decline)]
group_4_pfu_raw <- group_4[!group_4%in%c(exclude_pfu_raw, remove_raw_decline)]


# c(group_1_pfu_raw, group_2_pfu_raw, group_3_pfu_raw)
df_raw <- data.frame("id_sub"=character(), "win_length"=integer()) 
for (i in c(group_1_pfu_raw, group_2_pfu_raw)){
  pfu_i <- df_pfu$pfu[df_pfu$obs_id==i]
  # print(pfu_i[!is.na(pfu_i)])
  # first day
  pfu_first <- which(pfu_i>1)[1]
  pfu_last <- tail(which(pfu_i>1),1)
  df_raw <- rbind(df_raw, data.frame("id_sub"=unique(df_pfu$id_sub[df_pfu$obs_id==i]), "last_pfu"=pfu_last))
}

df_sym <- readRDS("IncidentSymptoms_v2.RDS")
pfu_raw <- merge(df_raw, df_sym, by="id_sub") %>%
  filter(include_sym) %>%
  filter(!is.na(as.double(contact_onset))) %>%
  mutate(recovery_from_sym = last_pfu - as.double(contact_onset))

# unvacc_from_sym <- pfu_raw
# vacc_from_sym <- pfu_raw



wilcox.test(unvacc_from_sym$recovery_from_sym,
            vacc_from_sym$recovery_from_sym)

quantile(pfu_raw$recovery_from_sym, prob=0.5)
quantile(pfu_raw$recovery_from_sym, prob=0.25)
quantile(pfu_raw$recovery_from_sym, prob=0.75)

quantile(vacc_from_sym$recovery_from_sym, prob=0.5, type=3)
quantile(vacc_from_sym$recovery_from_sym, prob=0.25, type=3)
quantile(vacc_from_sym$recovery_from_sym, prob=0.75, type=3)

quantile(unvacc_from_sym$recovery_from_sym, prob=0.5) 
quantile(unvacc_from_sym$recovery_from_sym, prob=0.25)
quantile(unvacc_from_sym$recovery_from_sym, prob=0.75)


# Vaccinated
# 5 (5, 6) N=14

# All 
# 6 (5, 7) N=27

# Unvaccinated
# 7 (5, 8) N = 13

# p = 0.14

# write.csv(df_raw$id_sub, "pfu_from_sym_all.csv", row.names = FALSE)
write.csv(vacc_from_sym$id_sub, "pfu_from_sym_vacc.csv", row.names = FALSE)
write.csv(unvacc_from_sym$id_sub, "pfu_from_sym_unvacc.csv", row.names = FALSE)



