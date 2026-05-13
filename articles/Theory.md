# Information-Based Network Inference

One promising area that is currently underutilized in mainstream
biology/bioinformatics, despite being computationally efficient and
solvable on a standard PC, is **Information-Based Network Inference
using “Part Mutual Information” (PMI)**.

While biologists often use correlation or standard Mutual Information
(MI) to build gene regulatory networks, these methods struggle with
“indirect dependencies” (e.g., if Gene A regulates Gene B, and Gene B
regulates Gene C, a standard tool might wrongly show a direct link
between A and C).

A common limitation mentioned in most recent bioinformatics publications
(such as those dealing with single-cell RNA sequencing or “multi-omics”
integration) is the **“curse of dimensionality”** and the
**false-positive rate** in network reconstruction.

Most current tools use “Partial Correlation,” which only works if the
relationship is linear (a straight line). However, biological systems
are non-linear. Researchers often admit their networks contain many
“false edges” because they cannot efficiently calculate non-linear
conditional independence on a standard computer.

**Part Mutual Information (PMI)** is a relatively new statistical metric
designed to identify the “direct” strength between two variables (like
two genes) while stripping away the influence of all other variables.

#### Why it’s better:

1.  **Non-linear:** Unlike standard correlation, it catches complex
    biological “curves.”
2.  **No “False Links”:** It specifically quantifies only the direct
    connection.
3.  **Low Power Requirements:** It does not require a supercomputer or a
    GPU. It uses basic math (logarithms and probability distributions)
    that a standard laptop CPU can calculate in seconds.

It is possible to use publicly available data from the **TCGA (The
Cancer Genome Atlas)** or **GEO (Gene Expression Omnibus)**, using a
*Step-by-Step Approach:*

1.  **Download Data:** Grab a simple CSV file of gene expression levels
    for a specific disease (e.g., Lung Cancer).

2.  **Define the Variables:** Let $`X`$ be Gene 1, $`Y`$ be Gene 2, and
    $`Z`$ be the “background” genes.

3.  **Calculate PMI:** Instead of the standard Mutual Information
    formula:

    ``` math
    I(X; Y) = \sum \sum P(x, y) \log \frac{P(x, y)}{P(x)P(y)}
    ```

    The PartMI package already implements the **PMI** formula, which
    uses a “partial” joint probability. This allows you to see if $`X`$
    and $`Y`$ are actually talking to each other, or if they just look
    like they are because $`Z`$ is talking to both of them.

### Why this hasn’t been widely implemented yet (I think!)

> Most biologists rely on “off-the-shelf” software like **WGCNA**
> (Weighted Gene Co-expression Network Analysis). WGCNA is old and uses
> simple correlation. So, creating a user-friendly, lightweight software
> tool that replaces WGCNA with **PMI** would be a major contribution to
> the field.

It solves the “limitation of indirect edges” frequently cited in papers,
it uses open data, and it runs perfectly on a normal home/office PC.
