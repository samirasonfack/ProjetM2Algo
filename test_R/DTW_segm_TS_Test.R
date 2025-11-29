seq_1 <- c(0, 1, 3, 1.2, 3.2, 1.4, 2, 1.3, 4,2, 0.5,0.5)
#n<- length(seq_1)
seq_2 <- c(0, 2, 4, 2.2, 3.2, 1.4, 0, 1.3, 3, 6, 1.5,0.5)
segment_size <- 4
step <- 2 
# We have segment_size > step, therefore we will have overlapping segments
# we chose this in order to avoid the information loss

segements_1 <- sliding_windows(seq_1, segment_size, step)
normalized_segments_1 <- normalize_segment(segements_1)

#print(n)
#print(segment_size)
#print(step)
#print(length(segements_1))

segements_2 <- sliding_windows(seq_2, segment_size, step)
normalized_segments_2 <- normalize_segment(segements_2)


plot(seq_1, type="l", col="black", lwd=2, ylim=range(c(seq_1, seq_2)),
    main="Séries 1 et 2 avec segments colorés", ylab="Valeur", xlab="Temps")

lines(seq_2, col="darkgrey", lwd=2)


cols_1 <- rainbow(length(normalized_segments_1))
for (i in seq_along(normalized_segments_1)) {
start <- (i-1)*step + 1
end <- start + segment_size - 1
lines(start:end, normalized_segments_1 [[i]], col=cols_1[i], lwd=2)
}

cols_2 <- rainbow(length(normalized_segments_2), start = 0.5, end = 1) 
for (i in seq_along(normalized_segments_2)) {
start <- (i-1)*step + 1
end <- start + segment_size - 1
lines(start:end, normalized_segments_2[[i]], col=cols_2[i], lwd=2, lty=2) 
}

# Dtw computation
dtw_results <- dtw_segment(normalized_segments_1, normalized_segments_2)
# Display DTW results
print(dtw_results)
