# Implemented ML architectures

These figures depict prediction flow in the implemented models, with the fitted components in purple. They are schematics, not new experimental results.

- [b16 PNG](b16_architecture.png) / [SVG](b16_architecture.svg)
- [Context model PNG](context_model_architecture.png) / [SVG](context_model_architecture.svg)
- [Two-page PDF](ml_architectures.pdf)

b16 uses a rank-128 signed ridge operator and global L1 amplitude matching. The operator acts on 17,742 measured gene coordinates; the complete output has 18,533 genes. Challenge and KOLF evaluation targets were excluded from operator training. Mean-transfer inputs can use the challenge targets' observed responses in source datasets. b17 changes the input to a weighted mean.

The context model uses standardized source-response components z and control components h. Standardized features [z, z×h] predict a low-rank correction that is decoded and added to a refitted b16-style denoiser. The displayed 48→16 regression is the rank-16, two-control-component example. Rank 8 gives 24→8; the actual ranks can be smaller if training data require it. The no-context ablation uses z alone. Source-only baseline predictions provide the new model's training residuals; all fitting and selection exclude the outer cell line.

Source implementations: `scripts/cross_source_residual.py`, `scripts/apply_cross_source_residual.py`, and `scripts/context_interaction.py`. Results: [completed experiment](../context_interaction_2026-09-15/RESULTS.md).
