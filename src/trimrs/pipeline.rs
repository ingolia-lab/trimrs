use serde::{Serialize, Deserialize};
use crate::qual_trim::{QualTrimEnd, QualTrimStart};
use crate::trimming::*;

#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Clone, Serialize, Deserialize)]
pub struct PipelineConf {
    qual_trim_start: QualTrimStartConf,
    qual_trim_end: QualTrimEndConf,
}

#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Clone, Serialize, Deserialize)]
pub struct QualTrimEndConf {
    threshold: u8,
}

impl QualTrimEndConf {
    pub fn qual_trim_end(&self) -> QualTrimEnd {
        QualTrimEnd::new(self.threshold)
    }
}

#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Clone, Serialize, Deserialize)]
pub struct QualTrimStartConf {
    threshold: u8,
}

impl QualTrimStartConf {
    pub fn qual_trim_start(&self) -> QualTrimStart {
        QualTrimStart::new(self.threshold)
    }
}

#[derive(Eq, PartialEq, Ord, PartialOrd, Hash, Clone, Serialize, Deserialize)]
pub struct Pipeline {
    qual_trim_start: Option<QualTrimStart>,
    qual_trim_end: Option<QualTrimEnd>,
}

impl Pipeline {
    pub fn process(&self, trimming: &mut Trimming) {
        if let Some(qts) = &self.qual_trim_start {
            qts.trim(trimming);
        }

        if let Some(qte) = &self.qual_trim_end {
            qte.trim(trimming);
        }
    }
}

