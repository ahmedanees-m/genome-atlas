"""Selection decision support — rank editors by use-case criteria."""
from __future__ import annotations

from dataclasses import dataclass

from genome_atlas.api import Atlas, EditorRecommendation


@dataclass
class UseCaseProfile:
    cell_type: str
    edit_type: str
    cargo_size_bp: int
    delivery: str
    prefer_dsb_free: bool


class SelectionEngine:
    def __init__(self, atlas: Atlas):
        self.atlas = atlas

    def rank(self, cell_type: str, edit_type: str, cargo_size_bp: int,
             delivery: str, prefer_dsb_free: bool,
             top_k: int = 5) -> list[EditorRecommendation]:
        profile = UseCaseProfile(cell_type, edit_type, cargo_size_bp,
                                 delivery, prefer_dsb_free)
        systems = self.atlas.systems()
        scored = [self._score_system(s, profile) for _, s in systems.iterrows()]
        scored = [s for s in scored if s is not None]
        scored.sort(key=lambda r: -r.pen_score)

        if not scored:
            return [EditorRecommendation(
                system="NO_MATCH",
                pen_score=0.0,
                aav_fit=False,
                mechanism="N/A",
                reasoning=["No system in ATLAS matches the constraints. "
                           "Consider relaxing cargo size or delivery requirements."]
            )]
        return scored[:top_k]

    def _score_system(self, sys_row, profile) -> EditorRecommendation | None:
        from genome_atlas.utils.size import system_total_size_aa

        reasoning = []
        name = sys_row.get("name", "")
        name_lower = name.lower()
        mech = sys_row.get("mechanism_bucket", "UNKNOWN")
        total_aa = system_total_size_aa(self.atlas, sys_row["node_id"])

        # ========== DSB AVOIDANCE (context-aware) ==========
        if profile.edit_type in ("deletion", "SNV") and profile.cargo_size_bp <= 10:
            # Simple edits: nucleases are established and efficient
            if mech == "DSB_NUCLEASE":
                s_dsb = 0.95
                reasoning.append("Nuclease: optimal for simple edits (deletion/SNV)")
            elif mech == "DSB_FREE_TRANSEST_RECOMBINASE":
                s_dsb = 0.7
                reasoning.append("DSB-free: capable but over-engineered for simple edits")
            else:
                s_dsb = 0.5
        elif profile.edit_type == "insertion" and profile.cargo_size_bp <= 200:
            # Small insertions: prime editors or compact nucleases + HDR
            if "prime" in name_lower:
                s_dsb = 1.0
                reasoning.append("Prime editor: purpose-built for small insertions")
            elif mech == "DSB_NUCLEASE":
                s_dsb = 0.85
                reasoning.append("Nuclease + HDR: well-established for small insertions")
            elif mech == "DSB_FREE_TRANSEST_RECOMBINASE":
                s_dsb = 0.75
                reasoning.append("DSB-free: capable alternative")
            else:
                s_dsb = 0.5
        elif profile.edit_type == "insertion" and profile.cargo_size_bp > 1000:
            # Large insertions: DSB-free mechanisms strongly preferred
            if mech == "DSB_FREE_TRANSEST_RECOMBINASE":
                s_dsb = 1.0
                reasoning.append("DSB-free: required for large cargo insertion")
            elif mech == "TRANSPOSASE":
                s_dsb = 0.9
                reasoning.append("Transposase: natural large-cargo mechanism")
            else:
                s_dsb = 0.2
                reasoning.append("Nuclease: poor fit for large cargo insertion")
        else:
            # Default / moderate cargo
            if mech == "DSB_FREE_TRANSEST_RECOMBINASE":
                s_dsb = 0.85
                reasoning.append("DSB-free mechanism")
            elif mech == "TRANSPOSASE":
                s_dsb = 0.4
                reasoning.append("Transposase: partial DSB via DDE")
            elif mech == "DSB_NUCLEASE":
                # C1: Cargo-tiered DSB penalty
                if profile.edit_type == "insertion" and profile.cargo_size_bp <= 1000:
                    s_dsb = 0.55
                    reasoning.append("DSB nuclease acceptable for small-insert HDR")
                elif profile.edit_type == "SNV":
                    s_dsb = 0.60
                    reasoning.append("DSB nuclease acceptable for precise SNV")
                else:
                    s_dsb = 0.20
                    reasoning.append("DSB-dependent; host repair required")
            else:
                s_dsb = 0.3

        # ========== DELIVERY / AAV FIT ==========
        if profile.delivery == "AAV":
            # C3: Unknown size guard
            if total_aa == 0:
                s_aav = 0.1
                reasoning.append("Size unknown; AAV compatibility uncertain")
            elif total_aa <= 600:
                s_aav = 1.0
                reasoning.append(f"Ultra-compact ({total_aa} aa); ideal for AAV")
            elif total_aa <= 900:
                s_aav = 0.9
                reasoning.append(f"Fits single AAV ({total_aa} aa <= 900)")
            elif total_aa <= 1200:
                s_aav = 0.7
                reasoning.append(f"Tight AAV fit ({total_aa} aa); may need optimization")
            elif total_aa <= 2000:
                s_aav = 0.5
                reasoning.append(f"Borderline AAV ({total_aa} aa); split-intein or dual AAV")
            else:
                s_aav = 0.1
                reasoning.append(f"Exceeds AAV capacity ({total_aa} aa > 2000)")
        else:
            s_aav = 0.9  # mRNA/LNP are cargo-agnostic

        # ========== CARGO SIZE FIT ==========
        # C2: SNV-specific prime editor boost (handled first, before edit-type branching)
        if profile.edit_type == "SNV" and "prime" in name_lower:
            s_cargo = 1.0
            reasoning.append("Prime editor: purpose-built for SNV correction")
        elif profile.edit_type == "insertion":
            if "prime" in name_lower:
                s_cargo = 1.0 if profile.cargo_size_bp <= 200 else 0.15
                if s_cargo < 0.5:
                    reasoning.append(f"Prime editing limited to ~200 bp (requested {profile.cargo_size_bp})")
                else:
                    reasoning.append(f"Prime editing: <= 200 bp supported")
            elif any(k in name_lower for k in ["recombinase", "bridge", "passige", "integrase"]):
                s_cargo = 1.0 if profile.cargo_size_bp <= 50_000 else 0.4
                if s_cargo < 0.5:
                    reasoning.append(f"Recombinase cargo may exceed safety limit ({profile.cargo_size_bp} bp)")
                else:
                    reasoning.append(f"Recombinase/integrase: large cargo OK")
            elif any(k in name_lower for k in ["cast", "transposon", "transposase"]):
                s_cargo = 1.0 if profile.cargo_size_bp <= 50_000 else 0.4
                reasoning.append(f"CAST/transposase: natural large-cargo capability")
            elif mech == "DSB_NUCLEASE":
                s_cargo = 1.0 if profile.cargo_size_bp <= 1000 else 0.2
                if s_cargo < 0.5:
                    reasoning.append(f"Nuclease+HDR: limited to short inserts ({profile.cargo_size_bp} bp)")
                else:
                    reasoning.append(f"Nuclease+HDR: short cargo OK")
            else:
                s_cargo = 0.5
        else:
            s_cargo = 0.9  # Deletions and SNVs are cargo-agnostic (unless prime editor caught above)

        # ========== CELL TYPE COMPATIBILITY ==========
        s_cell = 0.5
        human_cell_types = ("HEK293T", "K562", "Jurkat", "HSC", "HepG2", "iPSC")
        if profile.cell_type in human_cell_types:
            if any(k in name_lower for k in ["spcas9", "cas12a", "cas12f"]):
                s_cell = 0.95
                reasoning.append("Established human-cell editor (Cas9/Cas12)")
            elif "prime" in name_lower:
                s_cell = 0.9
                reasoning.append("Prime editing validated in human cells")
            elif any(k in name_lower for k in ["cast", "bridge", "evocast"]):
                s_cell = 0.85
                reasoning.append("Published activity in human cells")
            elif "fanzor" in name_lower:
                s_cell = 0.7
                reasoning.append("Fanzor: emerging human-cell data")
            elif any(k in name_lower for k in ["cre", "bxb1"]):
                s_cell = 0.8
                reasoning.append("Recombinase: well-established in mammalian cells")
            else:
                s_cell = 0.6

        # ========== COMPOSITE SCORE ==========
        # Adjust weights based on use case
        if profile.edit_type in ("deletion", "SNV") and profile.cargo_size_bp <= 10:
            # Simple edits: cell-type experience and AAV fit matter most
            weights = {"dsb": 0.2, "aav": 0.35, "cargo": 0.1, "cell": 0.35}
        elif profile.delivery == "AAV" and total_aa > 900:
            # AAV size-constrained: delivery score gets boosted
            weights = {"dsb": 0.25, "aav": 0.4, "cargo": 0.2, "cell": 0.15}
        else:
            weights = {"dsb": 0.3, "aav": 0.3, "cargo": 0.2, "cell": 0.2}

        if profile.prefer_dsb_free:
            weights["dsb"] = min(0.45, weights["dsb"] + 0.1)
            weights["aav"] = max(0.2, weights["aav"] - 0.05)
            weights["cargo"] = max(0.15, weights["cargo"] - 0.05)

        pen_score = (weights["dsb"] * s_dsb + weights["aav"] * s_aav +
                     weights["cargo"] * s_cargo + weights["cell"] * s_cell)

        return EditorRecommendation(
            system=name,
            pen_score=round(pen_score, 3),
            aav_fit=(s_aav >= 0.8),
            mechanism=mech,
            reasoning=reasoning,
        )
