no_folate_neural <- 631            # Neural tube defects (No folate)
folate_neural <- 24                     # Neural tube defects (Folate)
no_folate_premature <- 727    # Premature births (No folate)
folate_premature <- 563          # Premature births (Folate)

# Attributable Risk (AR) calculations
AR_neural <- no_folate_neural - folate_neural
AR_premature <- no_folate_premature - folate_premature

cat("Attributable Risk for Neural Tube Defects:", AR_neural, "per 100,000\n")
cat("Attributable Risk for Premature Births:", AR_premature, "per 100,000\n")
