# Cross-platform BANTER classification of false killer whale acoustic events

## Purpose

Working repository for 'Considerations when applying classification models across recording platforms: A case study with Hawaiian false killer whales' to be presented at DCLDE 2024 - Rotterdam, The Netherlands

## Links

[Cross-plaform BANTER Classification](https://drive.google.com/drive/u/0/folders/1AWEu-Z4wjEYvwvZVuYcyZXDN4tl2QFaP) (Google Drive Folder for large data files)
[Team notes](https://docs.google.com/document/d/1iQLqO2aHOUdpuAUZKB0HAanaTyXcnns8tLfGeC-rzr4/edit)
[ClassyPc repository](https://github.com/ybarkley/ClassyPc) (Related work by Y. Barkley)

## DCLDE Abstract

**Considerations when applying classification models across recording platforms: A case study with Hawaiian false killer whales**

*Selene Fregosi, Jenny Trickey, Megan Wood, Yvonne Barkley, Jennifer McCullough, Anne Simonis, Erin Oleson*

An effective two-step detection and classification approach was developed for identifying false killer whales (Pseudorca crassidens) from towed array recordings (DCLDE 2022 dataset). This approach uses Pamguard for detection of clicks, whistles, and burst pulses and then uses the banter R package for classification using a two-stage random forest approach. Events are classified as either ‘Pc’ (false killer whales) or ‘UO’ (unidentified odontocete). It would be useful to apply this approach to additional acoustic datasets that are not visually verified. But applying a model trained on one recording system will likely not result in robust classifications. Alternatively, a new model could be trained for each recording system, but in many cases would lack visual verification of species. We investigated the tradeoffs of applying a classification model trained on a different recording system, but with visually verified events, to a model trained on the same acoustic system but relying on expert assessment of manual species identification from recordings alone. We applied the model trained on towed array data (towed array model) to two additional datasets: a miniaturized HARP (High-frequency Acoustic Recording Package) deployed on longline fishing gear and a PMARXL (Passive Miniaturized Acoustic Recorder XL) installed in an underwater glider. Two trained analysts independently assessed over 200 cetacean events on both the longline and glider datasets to serve as the ‘ground truth’ dataset. Separate models were then trained for both the longline and glider datasets. Each model (towed array, longline, and glider) was subsequently applied to each dataset to predict species. Accuracy of the model predictions were assessed at three classification thresholds (80%, 60%, and 50%). Finally, a combined model was trained on all datasets and compared to the individual-recorder models. Preliminary results show that the accuracy varies by recorder type and is sensitive to detection and filtering steps implemented before model training as well as the number of signals of each type (clicks, whistles, and burst pulses) present in each event. Additional model improvement could likely be achieved through exploration of additional factors, such as a larger sample size or more diverse recording conditions. Because these automated approaches are invaluable in processing large amounts of passive acoustic data, understanding even small differences in model performance across recorder types is essential.
