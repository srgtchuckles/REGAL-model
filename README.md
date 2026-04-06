This is a parametric mixure cure model that aims to estimate the results of the Phase III REGAL trial. This trial compares Galinpepimut-S(GPS), a WT1 peptide vaccine, and Best Available Therapy(BAT) in Acute Myeloid Lukemia in Second Complete Remission(AML CR2). Currently, the BAT for AML CR2 would most likely be Venetoclax(VEN) and/or Azacitidine(AZA). In order for this trial to succeed, the hazard ratio(HR) must be <0.636 to indicate GPS incures a survival benefit of 36.4% compared to the BAT. 

Why do I use a parametric mixture cure model? 

Just using a Weibull or Cox model assumes that every single patient will eventually have an event, or the cumulative hazard will reach ∞ as time reaches ∞. However, trial data from the phase II study of GPS in AML CR2 shows that the Kaplan-Meier(KM) curve stops declining and flattens. This cannot happen in a standard survival model; it can't represent that. The cured(fraction π) people in this model have 0 hazard from this disease and won't die from AML CR2. The susceptible(fraction 1 - π) patients will eventually relapse and die from the disease. Their survival follows a standard parametric Weibull distribution. However, we don't know which group any patient is in. For example, a patient alive at month 60 could be cured, or just a long-tailed susceptible patient that hasn't relapsed yet. Because we can't observe group membership directly, I used the law of probability across both groups.

    S_weibull(t) = exp(-(t/λ)^k)
    S(t) = P(alive at t | cured) × P(cured)
        + P(alive at t | susceptible) × P(susceptible)

        = π  +  S_weibull(t) × (1 - π)

        = π + (1 - π) × S_weibull(t)

π = height of the long-run plateau

λ = the Weibull scale, found from the mOS of the non-cured patients

k = the Weibull shape, controls whether hazard is increasing(k>1), constant(k=1), or decreasing(k<1). In AML, we use k>1 because the hazard increases over time

Usually, for π, we can fit it statistically from patient-level data. However, because we don't have patient level data, we derived it from the biological characteristics and responses from the phase II population. 
