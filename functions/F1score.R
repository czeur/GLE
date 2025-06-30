# Function to calculate F1 score from a true Theta and an estimated Theta

F1score <- function(Theta, Theta_hat) {
    d <- dim(Theta)[1]
    if (d != dim(Theta_hat)[1]){
        stop("\n The two matrices do not match in dimensions!")
    }


      # Extract off-diagonal elements only
    Theta_off_diag <- Theta[upper.tri(Theta)]
    Theta_hat_off_diag <- Theta_hat[upper.tri(Theta_hat)]
  
    # Convert to binary indicators: 1 if there's an edge, 0 if not
    Theta_binary <- ifelse(Theta_off_diag != 0, 1, 0)
    Theta_hat_binary <- ifelse(Theta_hat_off_diag != 0, 1, 0)

    # Calculate true positives, false positives, and false negatives
    true_positives <- sum(Theta_binary == 1 & Theta_hat_binary == 1)
    false_positives <- sum(Theta_binary == 0 & Theta_hat_binary == 1)
    false_negatives <- sum(Theta_binary == 1 & Theta_hat_binary == 0)

    # Calculate precision and recall
    precision <- ifelse(true_positives + false_positives > 0, true_positives / (true_positives + false_positives), 0)
    recall <- ifelse(true_positives + false_negatives > 0, true_positives / (true_positives + false_negatives), 0)

    # Calculate F1 score
    F1 <- ifelse(precision + recall > 0, 2 * (precision * recall) / (precision + recall), 0)

    return(F1)
}


