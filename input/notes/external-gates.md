# External gates: VCC 2026 rules and the arXiv/bioRxiv route

Task 1 of `input/ideas/zero-shot-limits-context-complexity-plan.md`. Compiled 2026-11-05, the day of the VCC final-entry
deadline (23:59 UTC). Nothing has been sent, submitted or registered. Every quotation below is verbatim.

**How the VCC text was retrieved.** virtualcellchallenge.org renders its pages client-side, so WebFetch and Exa only
saw "Loading…". The text below was pulled from the site's own Next.js bundles on 2026-11-05: the rules page chunk
`/_next/static/chunks/2622-63c808ba4137a5e1.js`, the terms page, the FAQ page, and the layout banner chunk
`/_next/static/chunks/app/(public)/layout-44818b0882237a09.js`. It is the same text the browser renders at the URLs
cited. Before quoting any of it in the paper, open each URL in a browser and re-check it.

---

## 1. VCC 2026 rules on publication, data use and leaderboard data

### 1.1 Documents found

| Document | URL | Version |
|---|---|---|
| Virtual Cell Challenge 2026 Rules | https://virtualcellchallenge.org/rules | "As revised on September 16, 2026" |
| Website Terms of Use (part of the binding contract) | https://virtualcellchallenge.org/terms | "Effective Date: August 4, 2026" |
| FAQ | https://virtualcellchallenge.org/faq | undated |
| Arc announcement | https://arcinstitute.org/news/virtual-cell-challenge-2026 | 2026-08-20 |
| Cell commentary, "Virtual Cell Challenge 2026: Benchmarking zero-shot generalization across cellular contexts" | https://www.cell.com/cell/fulltext/S0092-8674(26)00931-1 | published 2026-09-17 |

The rules make the Terms binding: *"Please note these rules incorporate by reference an [Arbitration Agreement] … and the
Virtual Cell Challenge Website [Terms of Use] (“Terms of Use”). The Arbitration Agreement and Terms of Use form part of
the binding contract governing participation in the Challenge. If there is any conflict between these rules and any
other guidance about the Challenge, these rules control."* (https://virtualcellchallenge.org/rules, introduction)

Clarifications go out by email: *"If you have any questions about these rules, please email
help@virtualcellchallenge.org. Any clarifications to the rules will be announced to all entrants via email."* (same).
**Action for the user:** search your registration inbox for any clarification emails since 2026-08-04. I did not read
your mailbox.

### 1.2 Clauses on publishing methods or results

**FAQ: publishing (the only explicit statement).** https://virtualcellchallenge.org/faq, "Legal" section:
> "I'm a researcher or student. Can I publish my results?" — "We encourage publishing and sharing research with the
> scientific community, if that is allowed by your employer or academic institution."

> "Can I share my progress on social media?" — "Feel free to share about your experience."

The FAQ has no legal force where it conflicts with the rules, since the rules control. Here, though, the rules say
nothing on the subject, so the two do not conflict.

**Rules: publication by Arc of finalists' methods.** https://virtualcellchallenge.org/rules#finalists-winners-prizes:
> "If your team is identified as a finalist, in order to be eligible for a prize, you … must: To help others
> participate in advancing the field, submit a high-level description, which we will publish, of how you created your
> winning entry. You need not disclose code or model weights and the description does not need to be detailed enough
> for others to be able to directly reproduce your approach, although we encourage as much transparency as you are
> comfortable providing."

This is a disclosure duty on finalists. It is not an embargo, and nothing makes prize eligibility depend on keeping
methods unpublished.

**Rules: collusion during the competition.** https://virtualcellchallenge.org/rules#responsibility-to-others:
> "You may not … collude with other teams to share models, results or strategies in a manner that undermines the
> integrity of the Challenge; attempt to reverse engineer, tamper with or otherwise gain unauthorized access to another
> team's work; or engage in any conduct that is detrimental or unfair to us or any other team, as we may determine in
> our sole discretion."

This concerns conduct during the competition. A preprint posted after the 2026-11-05 final-entry deadline cannot
affect any entry.

**Timeline and final results.**
- Rules, #submitting-final-entries: *"You must upload the final entry to the Challenge website by 11:59pm UTC on
  November 5, 2026."* … *"There will be no leaderboard displayed for final entries from October 22, 2026 to November 5,
  2026."* … *"A final leaderboard, based only on each team's last final entry, will be released when the winners are
  announced."*
- Site banner (layout chunk): *"The 2026 Virtual Cell Challenge is now closed. We will announce winners and publish the
  final leaderboard in late November."*
- FAQ: *"Winners will be announced in late November."*
- Arc news page: *"Mid-late November: Winners announced"*.

**Embargo.** No clause in the rules, the Terms, the FAQ, the Arc announcement or the Cell commentary forbids or delays
publication before the final results are announced. None ties prize eligibility to non-publication either. I searched
the whole rules text for publish, embargo, preprint, paper, confidential and announce.

### 1.3 Clauses on the official data (validation and final control cells)

**Rules.** The rules describe the data but grant no licence and set no restriction on publishing with it:
> "At the start of the Challenge we will release the non-targeting control (NTC) profiles for the three validation cell
> lines, and the validation perturbations (300 per cell line)." (#submitting-leaderboard-entries)
> "The final-test perturbations and the non-targeting control (NTC) profiles for the three held-out test cell lines will
> be released on October 22, 2026." (#submitting-final-entries)

**Terms of Use: "Content", which covers downloaded data.** https://virtualcellchallenge.org/terms:
> "Our Services may include text, graphics, data, instructions, articles, photos, images, illustrations and other users'
> content (collectively, the “Content”). The Content is protected by copyright and other intellectual property laws, and
> forms part of our Services, even if you download it. Your use of our Content is subject to these Terms, including the
> rules and restrictions in the following section."

Under "Rules and Restrictions": *"You agree that you will not do any of the following:"*
> "Use your account for any purpose other than participating in the Challenge;"
> "Use, copy, reproduce, modify, translate, publish, broadcast, transmit, distribute, perform, upload, display, license,
> sell, create derivative works based on, or otherwise exploit for any purpose any of our Services except as expressly
> permitted by these Terms;"
> "Remove, alter or obscure any copyright, trademark, service mark or other proprietary rights notices …"

> "Arc Research Institute and Virtual Cell Challenge are trademarks and/or registered trademarks of Arc Research
> Institute. You may not use these trademarks without Arc Research Institute's prior written permission …"

**No data licence found.** The downloaded validation bundle
(`vcc2026_data/raw/vcc_official_controls_refetch/extracted/manifest.json`, panel `vcc2026-val-1`) contains no licence
file. The `vcc` CLI is MIT-licensed, but that covers code only
(`vcc_cli-0.3.0.dist-info/licenses/LICENSE`, "Copyright (c) 2026 Arc Research Institute"). The 2025 dataset was later
released through the Arc Virtual Cell Atlas (*"The full 2025 dataset … is now downloadable through the Arc Virtual
Cell Atlas"*, https://arcinstitute.org/news/virtual-cell-challenge-2026). No release or licence has been announced
for the 2026 data.

**Training-data FAQ (about other data, not the official data):** *"As long as you have the appropriate permissions to
use the data, you may use any data to improve the performance of your predictions."*
(https://virtualcellchallenge.org/faq)

### 1.4 Clauses on leaderboard data and other teams

**Terms of Use: licence from entrants to Arc.** https://virtualcellchallenge.org/terms, "Participating in the
leaderboard and similar Services":
> "If you voluntarily choose to participate in any publicly accessible or participant-facing components of the Services,
> including any leaderboard … you grant us a worldwide, non-exclusive, royalty-free, sublicensable and transferable
> license to use, host, store, reproduce, modify, publish, display, perform and distribute any information, materials or
> content that you submit … for the purposes of operating and administering the Services and the Challenge. You
> acknowledge and agree that such content may be publicly visible."

Leaderboard rows are "other users' content" and therefore "Content" (definition quoted in §1.3). They fall under the
same "publish … create derivative works … except as expressly permitted" restriction. Nothing expressly permits
entrants to republish the leaderboard, and nothing forbids it either. The leaderboard page currently shows: *"Live
rankings are now frozen and will no longer update. Final standings will be shown in the Final Results tab after the
challenge ends."* (https://virtualcellchallenge.org/leaderboard)

**Mentioning other teams:** neither the rules nor the Terms say anything about this. The Terms' ban on "defamatory"
material applies only to material posted *to the Services*.

### 1.5 Verdict for a preprint posted on 2026-11-28

| Question | Verdict | Basis |
|---|---|---|
| (a) Report own validation-phase scores and rank (129 / 1,635) | **Allowed.** | FAQ: "We encourage publishing and sharing research … if that is allowed by your employer or academic institution." No embargo anywhere. The final-entry deadline has passed (2026-11-05), so the collusion clause cannot apply. Your scores are your own results. One condition: if you registered as a *university representative* (rules #participating-as-employee), check that UCSD / your institution allows it. |
| (b) Use the official validation and final NTC control cells | **Ambiguous. Ask the organisers.** | The rules are silent and there is no data licence. A literal reading of the Terms ("Content … even if you download it"; no "use … publish … create derivative works … except as expressly permitted"; account use only "for participating in the Challenge") could cover a paper's analyses. The FAQ's encouragement to publish "results" points the other way. Low-risk practice until Arc answers: **do not redistribute** any cells or count matrices (no Zenodo, no supplement), report only derived summary statistics, cite the Cell commentary, and point readers to virtualcellchallenge.org for access. Analysing the *final-phase* controls (D/E/F) before the winners are announced is not prohibited either, but ask about it in the same email. |
| (c) Aggregate public leaderboard data without team names | **Ambiguous, but low risk.** | No clause grants or denies it. Leaderboard rows are "Content" under the Terms. Aggregate statistics (score distribution, quantiles, the number of teams beating the context-mean baseline) with no team names involve no personal data and little copyright exposure. Still, include it in the question to Arc. If Arc says no, fall back to reporting only your own rank and the published team count. |
| Must publication wait until winners are announced (late November)? | **No rule requires it.** | No embargo clause exists. As a courtesy, and to be able to report final-phase scores, aim to post *after* the winners announcement. 11-28 to 12-01 probably satisfies this, but confirm on the day. |
| Prize eligibility | **Not affected.** | No clause links publication to eligibility. Finalists must in any case submit a description that Arc publishes. |
| Trademarks | **Take care.** | Use "Virtual Cell Challenge" only descriptively and with a citation. Do not use logos and do not imply Arc endorses the paper. |

### 1.6 Draft clarification email to the organisers (text only, not sent)

To: help@virtualcellchallenge.org (the contact address given in the rules)
Subject: Publication question, 2026 entrant (team "[TEAM NAME]")

> Hello VCC team,
>
> Thank you for running the 2026 Challenge. I am preparing a preprint (target late November / early December 2026,
> after the winners are announced) on the limits of zero-shot perturbation-response prediction in unseen cell lines,
> drawing on my participation. Before posting, I would like to confirm three points:
>
> 1. May I report my own validation-phase scores and leaderboard rank, and, once released, my final-phase score?
> 2. May I report analyses computed on the released NTC control profiles (validation contexts A–C and final contexts
>    D–F), such as summary statistics and model evaluations, without redistributing any cells or count matrices,
>    citing your Cell commentary and pointing readers to virtualcellchallenge.org for the data?
> 3. May I show aggregate statistics of the public validation leaderboard (e.g. the score distribution), with no team
>    names?
>
> If any of these needs adjusting, or if you would prefer that I wait for a specific date, I am happy to follow your
> guidance.
>
> Best regards,
> [Your Name]

---

## 2. arXiv endorsement for a first-time cs.LG submitter

### 2.1 Current policy (checked 2026-11-05)

- **Endorsement is required:** *"arXiv requires that users be endorsed before submitting their first paper to arXiv or
  a new category."* (https://info.arxiv.org/help/endorsement.html)
- **Policy change of 2026-01-21:** an institutional email is no longer enough on its own. A new submitter is now
  endorsed automatically only if they have **both** an academic/research email **and** prior authorship of an accepted
  paper in the relevant endorsement domain. Everyone else needs personal endorsement from an established author. arXiv
  says "institutional email addresses are no longer a sufficient credential"; authors already endorsed in a category
  stay endorsed. (https://blog.arxiv.org/2026/01/21/attention-authors-updated-endorsement-policy). This extends the
  Mathematics pilot of 2025-12-10 (https://blog.arxiv.org/arxiv/2025/12/10/updated-endorsement-policy-for-arxiv-mathematics).
  *What this means here:* if you have no prior cs paper submitted from your account, expect to need a personal
  endorser even with a ucsd.edu address.
- **Who can endorse:** *"Endorsers must have authored a certain number of papers within the endorsement domain of a
  subject area"*. Only papers *"submitted between three months and five years ago"* count, and the endorser must hold
  an active endorsement in that area. The help page gives no number for cs. The page says *"most high-level subject
  areas (e.g., hep-th, cond-mat, q-bio) are currently endorsement domains"* but does not say outright that all of cs is
  one domain. (https://info.arxiv.org/help/endorsement.html)
- **Etiquette:** *"it is inappropriate to email large numbers of potential endorsers at once"*; do not re-email the same
  person repeatedly. Endorsers are told *"If you feel uncomfortable about endorsing an author for any reason, do not do
  it."* (same page). **Contact one candidate at a time.**
- **Cross-list to q-bio.GN:** neither the endorsement page nor https://info.arxiv.org/help/cross.html says whether a
  cross-list needs separate endorsement. The cross-list page warns: *"It is rarely appropriate to add more than one or
  two cross-lists"* and *"Bad cross-lists will be removed."* A single q-bio.GN cross-list is consistent with that. If
  the submission form asks for q-bio endorsement, drop the cross-list rather than delay.
- **Related 2025–26 CS rules:**
  - From 2025-10-31, CS review articles and position papers need documented prior peer review
    (https://blog.arxiv.org/?p=1740, "Attention Authors: Updated Practice for Review Articles and Position Papers in
    arXiv CS Category"). The preprint must therefore read as a **research paper with new experiments**, not a
    position paper. Avoid "Position:" framing in the title and abstract.
  - From 2026-10-01, a rate limit applies: *"up to two submissions per calendar month"* and *"a limit of three total
    active submissions at any given time"*. AI use is allowed *"as long as it is disclosed"*
    (https://blog.arxiv.org/2026/10/01/updated-rate-limit-policy/). Disclose AI assistance in the acknowledgements.

### 2.2 How the endorsement code is obtained

1. Create or log in to an arXiv account (https://arxiv.org/user/), using the institutional email if you have one.
2. Start a new submission and choose **cs.LG** as the primary category. Because the account is not endorsed for cs,
   arXiv stops the submission and sends an **endorsement request email**. This email contains a **six-character
   alphanumeric endorsement code** and a link you can forward to prospective endorsers
   (https://info.arxiv.org/help/endorsement.html: "Check your email for an endorsement request email. This will include
   a link you can provide to prospective endorsers.").
3. Find a qualified endorser. On the arXiv abstract page of a related recent paper, click **"Which authors of this
   paper are endorsers?"** (requires login: `https://arxiv.org/auth/show-endorsers/<arXiv-id>`). This confirms
   eligibility before you write. The endorser's address is under "Submission history" on the abstract page. Prefer
   the public contact page listed below.
4. Send the code to **one** endorser. The endorser enters it at https://arxiv.org/auth/endorse. One positive
   endorsement per endorsement domain is enough.
5. Resume the submission. Plan gate: **if there is no endorsement by 2026-11-21, switch to bioRxiv and log the switch.**

### 2.3 Candidate endorsers (contact in this order, one at a time)

Before writing to anyone, confirm eligibility with the "Which authors of this paper are endorsers?" link. I could not
check it, because it needs login and the arXiv API returned HTTP 429.

| # | Name, affiliation | Related paper | Why they fit | Public contact page |
|---|---|---|---|---|
| 1 | **Han Liu**, Orrington Lunt Professor of Computer Science (and Statistics), Northwestern University; directs the MAGICS Lab and the Center for Foundation Models & Generative AI | Dibaeinia, …, **Han Liu**, Perera, Khan, *"Position: Virtual Cells Need Context, Not Just Scale"*, ICML 2026 (https://icml.cc/virtual/2026/poster/67166; preprint https://www.biorxiv.org/content/10.64898/2026.02.04.703804v1). Also *Cell-JEPA* (arXiv 2602.02093, https://arxiv.org/abs/2602.02093) | The position paper argues that context coverage, not scale, limits virtual cells. Your paper is the empirical, context-count side of that argument. He is a long-standing ML author with many NeurIPS/ICML/ICLR papers (lab page: "48 NeurIPS · ICML · ICLR papers, 2023–2026"), so he is very likely a cs endorser. Note: the position paper itself is on bioRxiv, not arXiv. | https://magics.cs.northwestern.edu/people.html ; https://www.mccormick.northwestern.edu/research-faculty/directory/profiles/liu-han.html |
| 2 | **Sungsoo Ahn**, Assistant Professor, KAIST Kim Jaechul Graduate School of AI (SPML Lab) | Piao, …, Park, **Ahn**, *"Learning Adaptive Perturbation-Conditioned Contexts for Robust Transcriptional Response Prediction"* (AdaPert), ICML 2026, PMLR 306 (https://arxiv.org/abs/2602.18885; https://proceedings.mlr.press/v306/piao26a.html). Also *"Progressive Multi-Agent Reasoning for Biological Perturbation Prediction"* (arXiv 2602.07408) | He is the senior author on a recent ML method for perturbation-response prediction. The arXiv author search lists 77 papers, with 6 of his 15 most recent being cs.LG primary (https://arxiv.org/search/?query=Sungsoo+Ahn&searchtype=author), so he is almost certainly a cs endorser. AdaPert itself is cs.CE. | https://sungsoo-ahn.github.io/ ; https://gsai.kaist.ac.kr/people/ |
| 3 | **Kun Zhang**, Professor, Carnegie Mellon University (Philosophy; affiliate, Machine Learning Dept.); also leads CLeaR at MBZUAI | Wu, Wershof, Schmon, Nassar, Osiński, Eksi, Yan, Stark, **Kun Zhang**, Graepel, *"PerturBench: Benchmarking Machine Learning Models for Cellular Perturbation Analysis"*, arXiv 2408.10609, **cs.LG primary, q-bio.GN + stat.ML cross-lists** (https://arxiv.org/abs/2408.10609; NeurIPS 2025 slides: https://neurips.cc/media/neurips-2025/Slides/121669.pdf) | PerturBench is a benchmark paper with the same category combination you plan to use, and it finds that "simpler architectures are generally competitive". His core area is transfer learning across heterogeneous environments and causal ML (CMU bio: "Domain adaptation/transfer learning; Learning in nonstationary/heterogeneous environments"), which matches the context sample-complexity framing. He is a highly established cs.LG author. | https://www.cmu.edu/dietrich/philosophy/people/faculty/zhang.html ; https://mbzuai-clear.github.io/ |

Reserves, if the three above decline or do not answer:
- **Aly A. Khan** (U Chicago; corresponding author of the position paper). The arXiv author search shows 4 cs papers in
  the 3-month-to-5-year window: 2112.13737, 2311.05659 and 2510.18122 (all cs.LG), and 2602.02093 (cs.CE).
  https://arxiv.org/search/?query=%22Aly+A.+Khan%22&searchtype=author
- **Yan Wu** (Altos Labs; first author and submitter of PerturBench).

Caveat: Altos Labs won a 2025 VCC prize and may have competed in 2026. To avoid any appearance of conflict, leave Yan
Wu until after the winners are announced.

The plan's overlap paper, Molina & Zhang, was deliberately not listed. Its authors are closest to your topic, but
the overlap gate in Task 2 is still open.

### 2.4 Draft endorsement-request email (≤150 words; two placeholders)

Send it to one candidate's public contact address. Attach the abstract and, if ready, the PDF draft.

> Subject: arXiv cs.LG endorsement request, zero-shot perturbation-response prediction
>
> Dear Professor,
>
> I am an independent researcher preparing my first arXiv submission and would be grateful for an endorsement in
> cs.LG. Your work on perturbation-response prediction is closely related, which is why I am writing to you
> specifically.
>
> The paper studies the limits of zero-shot perturbation-response prediction in unseen cell lines. It decomposes
> responses into transferable and line-specific parts, compares learned models with simple baselines against noise
> ceilings, and measures how performance scales with the number of source cell lines. It builds on my participation
> in the Arc Virtual Cell Challenge 2026 (validation rank 129 of 1,635). The abstract and draft are attached.
>
> If you are willing, my endorsement code is [ENDORSEMENT CODE]; it can be entered at https://arxiv.org/auth/endorse.
> I fully understand if you prefer not to.
>
> Thank you for your time.
>
> Best regards,
> [Your Name]

If you registered for VCC as a UCSD representative, replace "independent researcher" with your affiliation. Check that
the body stays at or under 150 words after any edit.

---

## 3. bioRxiv fallback (if there is no endorsement by 2026-11-21)

- **No endorsement is needed.** *"Authors wishing to deposit manuscripts must first register on the site and supply an
  affiliation so there is an entity providing oversight …"*; *"There is no charge for registration or article
  deposition."* Submit at https://submit.biorxiv.org (https://www.biorxiv.org/submit-a-manuscript). An affiliation is
  required, so use UCSD or another real affiliation.
- **Scope:** *"bioRxiv welcomes computational papers. These should be full research papers and include methodological
  details and results. Simple automated/computational analyses of public data … are generally not sufficient"*.
  Hypotheses without new data, opinion pieces and narrative reviews are out of scope
  (https://www.biorxiv.org/about/FAQ). As with arXiv, the paper must present new analyses, not a position.
- **Category:** *"Only one subject area can be selected for a preprint on bioRxiv."* Recommended: **Bioinformatics**.
  Alternatives are Genomics and Systems Biology. Article type: **New Results**
  (https://www.biorxiv.org/submit-a-manuscript).
- **Prior arXiv posting:** *"We strongly recommend that a preprint be posted on only one server … Manuscripts already
  present on other preprint servers will not be posted."* (https://www.biorxiv.org/about/FAQ, "Can I submit a manuscript
  to bioRxiv and another preprint server?"). arXiv and bioRxiv are **mutually exclusive** for this paper. Do not post
  to arXiv once the bioRxiv submission is made, and if it is already on arXiv, bioRxiv will refuse it.
- **Timing:** *"Preprints usually appear on bioRxiv within 72 hours; however, screening may take longer over a weekend
  or holiday"* (FAQ). US Thanksgiving is 2026-11-26, so to post by 12-01, submit by about 11-25.
- **Format and licence:** submit a single PDF (LaTeX must be converted to PDF). Licences offered: CC BY, CC BY-NC,
  CC BY-ND, CC BY-NC-ND, CC0, or no reuse. Posted preprints cannot be removed, only withdrawn. Describe any generative-AI
  use in the manuscript (FAQ).

---

## 4. Open actions for the user

1. Check your inbox for any VCC rule-clarification emails since 2026-08-04, then decide whether to send the §1.6
   question. Sending it this week leaves about three weeks for a reply.
2. Check which capacity you registered in for VCC (personal or university representative). If you registered as a
   university representative, confirm that UCSD allows the publication.
3. Create or log in to your arXiv account and start the cs.LG submission to trigger the endorsement code. Verify
   candidate 1 with "Which authors are endorsers?", then send the §2.4 email to that one person. Move to the next
   candidate only after about 5–7 days without a reply. Hard switch to bioRxiv on 11-21.
4. Whatever Arc answers on (b), do not put VCC cells or count matrices into any public supplement.
