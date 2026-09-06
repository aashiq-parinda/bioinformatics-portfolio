# Ashraf Khan

*Senior Full-Stack & Systems Engineer | Applied GenAI: LLM Fine-Tuning, Alignment & Inference Optimization*

Mumbai, Maharashtra, India | +91-8779559898 | [ashrafk.salim1@gmail.com](mailto:ashrafk.salim1@gmail.com) | [GitHub](https://github.com/aashiq-parinda) | [Hugging Face](https://huggingface.co/ashrafksalim) | [LinkedIn](https://linkedin.com/in/ashrafksalim)

---

## Experience

**Logicloop** — Jan 2026 – Present
*Senior Software & Systems Engineer — Mumbai, India*

- **[Lead Yatra](https://apps.apple.com/in/app/lead-yatra/id6758609992)** (Reliance IndusInd) [[Android](https://play.google.com/store/apps/details?id=com.rnlic.leadyatra&hl=en_IN)]: lead architect and engineer for a cross-platform Flutter application optimizing field sales lead capture and venue booking at enterprise scale — handling sensitive FinTech and banking data for IndusInd Bank's field operations, powered by GCP Cloud Run microservices and Google Maps Platform.
- Built the [MAXX](https://apps.apple.com/us/app/maxx/id6747980186) AI Voice Assistant with custom wake-word activation and multi-advisor conversational flows via Millis.ai, applying LLM-driven conversation design to a production financial services product.
- Led a team of 3 developers; enforced CI/CD automation and Level 2 VAPT security standards across enterprise financial flows, cutting deployment overhead by 40%.

*Also: Software Engineer at Logicloop (May 2024 – Oct 2024) — built the Flutter frontend for [HDFC Life InstaQuote](https://apps.apple.com/in/app/instaquote-2-0-by-hdfc-life/id6639613496), integrating insurance-pricing ML backend APIs.*

**Fitwell Technologies Inc.** — Jun 2025 – Jan 2026
*Full-Stack Developer & Data Custodian — Remote (VA, USA)*

- Integrated LLM-powered features and Amazon Connect voice AI workflows into the [Monitor Health AI](https://monitorhealth.ai/) platform on HIPAA-compliant AWS pipelines.
- Deployed a self-hosted Supabase/PostgreSQL backend on EC2 for PHI data isolation; maintained zero data breach posture as data custodian.

**QuantGen Private Limited** — Nov 2024 – Jun 2025
*Lead Software Engineer — India*

- Delivered an AI-integrated health platform ([FitWell AI](https://apps.apple.com/in/app/fitwell-all-in-one/id1641657795)) and an enterprise 24/7 operations platform ([Kokomo 24/7](https://apps.apple.com/us/app/kokomo24-7/id1504839321)) with AWS ECS/App Runner orchestration.

**Additional Experience**

- SoftProdigy, Remote USA — Software Engineer (Nov 2023 – Apr 2024): loan lending application workflows; gRPC client API libraries.
- AKcess Labs Pvt. Ltd., Remote UK — Senior Software Development Engineer (Jan 2023 – Nov 2023): Flutter + Flask hyperlocal commerce platform.
- We3.Tech, Maharashtra — Software Developer (Aug 2021 – Jan 2023): eKYC, Digio eSign, and banking API modules for Motilal Oswal MO Investor (SEBI/RBI-compliant).

---

## Technical Skills

- **LLM Optimization & Inference:** 4-bit NF4 Quantization, bitsandbytes, GGUF, Prefix KV-Cache / Prompt Caching, Torch float16, VRAM Profiling, T4/L4 GPU Deployment
- **Fine-Tuning & Alignment:** LoRA, QLoRA, DPO (Direct Preference Optimization), Bradley-Terry Reward Modeling, PEFT Adapters, DeBERTa-v3, Qwen2.5, Hugging Face transformers & trl, LEDGAR Corpus
- **Agentic Pipelines & RAG:** Multi-Stage DAGs, ChromaDB, BM25 + Dense Vector Hybrid Search, RRF, Cross-Model Context Handoff, LangChain, Prompt Engineering
- **Safety & Evaluation:** Multi-Tier Guardrail Catalogs, Adversarial/Jailbreak Defense, LLM-as-Judge (G-Eval), Structured Output (FSM Grammar Sampling & Logit Masking), p50/p90/p95/p99 Latency Benchmarking, F1/Precision/Recall, Pytest CI/CD
- **APIs & Serving:** FastAPI (Async), Server-Sent Events, REST, gRPC, Protobuf, WebSockets
- **ML Frameworks & Cloud:** PyTorch, TensorFlow, scikit-learn, OpenCV, AWS (EC2, S3, ECS, Fargate, App Runner), GCP, Docker, GitHub Actions CI/CD
- **Programming & Data:** Python (Expert), SQL/PostgreSQL, Flutter/Dart, React Native, Node.js

---

## Featured GenAI / ML Projects — [Full Portfolio](https://github.com/aashiq-parinda/genai-systems-portfolio)

### Enterprise GenAI Platform & Capacity Architecture [[Code]](https://github.com/aashiq-parinda/genai-systems-portfolio/tree/main/portfolio/Enterprise-GenAI-Platform-Capacity-Architecture)
Multi-tenant `/v1/chat` control plane for a Tier-1 enterprise (10M users, 100K DAU) with dynamic SLM/Frontier model routing, a zero-trust tool plane, and first-principles GPU capacity modeling from 10K to 1M concurrency. Engineered a 3-layer context caching hierarchy — Prefix KV-Cache (SHA-256 prompt hashing) and SLM→Frontier context-handoff compression — modeling 90%+ prefix token savings.
**Verified: 113.4x ROI, ₹70–90 Cr/yr infra cost avoidance, ₹120–150 Cr platform revenue modeled.**

### Customer Support Automation — Quantized LLM + Agentic Workflow [[Code]](https://github.com/aashiq-parinda/Customer-Support-Automation-Quantized-LLM-Agentic-Workflow) [[Live Demo]](https://huggingface.co/spaces/ashrafksalim/customer-support-quantized-llm-agent) [[Notebook]](https://colab.research.google.com/drive/1RIf9_bZAoqmB9rq6nNWciNmaSPZJSpOx?usp=sharing)
4-bit NF4-quantized 7B-class LLM agent with multi-tier safety guardrails and hybrid BM25 + dense-vector RAG, served via FastAPI SSE.
**Verified: 96.67% cost reduction, 65% VRAM reduction, 7.68x latency speedup, 0.94 guardrail F1, 93.75% accuracy preserved.**

### Contract Risk Review — LoRA Fine-Tuning & Multi-Agent Legal Reasoning [[Code]](https://github.com/aashiq-parinda/LLM-Fine-Tuning-Agentic-Reasoning-for-Contract-Risk-Review)
LoRA/PEFT-fine-tuned DeBERTa-v3/Qwen2.5 clause classifier on the LEDGAR corpus, with ChromaDB precedent-retrieval RAG and an agentic DAG routing critical clauses to human review.
**Verified: 1.000 macro F1, 0.12ms mean inference latency.**

### GenAI Core Algorithms, Evaluation & Alignment Masterclass [[Code / Notebooks]](https://github.com/aashiq-parinda/genai-systems-portfolio)
First-principles implementations in pure Python/NumPy: DPO closed-form loss & Bradley-Terry reward models, G-Eval LLM-as-Judge, FSM-based constrained decoding, BPE tokenizer, GQA with KV-Cache, RoPE, and a 58K-parameter Nano-Transformer trained in under 15s on CPU.

### Quantum Hardware Validation & Benchmarking Suite [[Code]](https://github.com/aashiq-parinda/Quantum-Hardware-Validation-Benchmarking-Suite)
Hardware-agnostic NISQ benchmarking engine verifying 6 core quantum-mechanical axioms with density-matrix noise simulation and Zero-Noise Extrapolation.
**Verified: all 6 axioms passed, +6.25% error recovery via ZNE.**

---

## Certifications

AI & Quantum Computing Mastery · Machine Learning with Python · Complete Quantum Computing Course · Healthcare API Compliance & HIPAA PHI Data Governance · FinTech Regulatory Standards & PCI-DSS Payment Processing
