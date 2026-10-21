# Data sources

All web sources were accessed on 10 October 2026. The same references, with full details, are in the
bibliography of the Supporting Information (SI).

## Calibration of the case study (base year 2024)

The numbers are in `model/calibration_2024.py`. Each quantity is taken from one official source; press
reports are used only where no official record exists.

| Used for | Source |
|---|---|
| Expressway (CT04) and QL5 toll revenue and traffic, 2024; number of BOT projects and stations | Vietnam Road Administration (Cuc Duong bo Viet Nam), "Cong khai so lieu thu phi va luu luong cac du an BOT" [Monthly toll revenue and traffic of the BOT road projects under the Ministry of Transport, January to December 2024; spreadsheet annex "Tong hop BOT 2025"], published with the notice for March 2025. <https://drvn.gov.vn/tin-tuc/hoat-dong-nganh-duong-bo/cuc-duong-bo-viet-nam-cong-khai-so-lieu-thu-phi-va-luu-luong-cac-du-an-bot-thang-3-2025.html> (accessed 10 October 2026) |
| Interest paid to VDB in 2024, debt to VDB, VDB share of equity, toll revenue excluding VAT | VIDIFI (Vietnam Infrastructure Development and Finance Investment Joint Stock Company), Consolidated financial statements for the year ended 31 December 2024, audited by NTV Auditing Company, 26 February 2025: income statement (interest expense), note V.13 (loans from VDB), note V.19 (owners' equity), note VI.1 (toll revenue). <https://vidifi.com.vn/upload/bctc/bctc-hop-nhat-vidifi-2024.pdf> (accessed 10 October 2026) |
| Full-trip expressway fares and QL5 fares per station (five vehicle classes) | VIDIFI, "Thong tin dieu chinh gia ve dich vu su dung duong bo" [Official tariff from 29 December 2023: expressway fares by section, including the whole tolled length, and single-trip fares per QL5 station, five vehicle classes], 21 December 2023. <https://evidifi.vn/vi/thong-tin-dieu-chinh-gia-ve-dich-vu-su-dung-duong-bo.html> (accessed 10 October 2026) |
| Vehicles per day through the QL5 stations and on the road (proxy phi_emp) | Bao Dau thau, "Nghich ly thu phi tai Quoc lo 5: Bo GTVT cung rat dau dau" [Operator's count: about 30,000 vehicles per day through the QL5 stations, 40,000 to 45,000 on the road], 17 October 2016. <https://baodauthau.vn/nghich-ly-thu-phi-tai-quoc-lo-5-bo-gtvt-cung-rat-dau-dau-post18447.html> (accessed 10 October 2026) |
| Share of revenue retained under manual collection (6%) | Tien Phong, "Thu phi tu dong khong dung: ai dang tao them cac rao can?" [ETC service share ca. 5% vs 6% retained under manual collection], 2019. <https://tienphong.vn/thu-phi-tu-dong-khong-dung-ai-dang-tao-them-cac-rao-can-post1091649.tpo> (accessed 10 October 2026) |
| National range of toll collection costs (3% to 12% of revenue) | Nha dau tu, "Tram thu phi khong dung: co tang thoi gian thu phi BOT?" [Toll-organisation cost 3-12% of revenue], 2017. <https://nhadautu.vn/tram-thu-phi-khong-dung-co-tang-thoi-gian-thu-phi-bot-d4751.html> (accessed 10 October 2026) |
| Fees for electronic toll collection (2.15% and 4.5% of revenue) | Tai chinh Doanh nghiep, "Yeu cau dung thu phi 3 tram BOT: ap dat vo ly va gay bat on" [Cam Thinh ETC contract: 2.15% vs the mandated 4.5% of revenue], 2019. <https://taichinhdoanhnghiep.net.vn/yeu-cau-dung-thu-phi-3-tram-bot-ap-dat-vo-ly-va-gay-bat-on-d8178.html> (accessed 10 October 2026) |
| The corridor went to electronic collection only on 1 June 2022 | Bao Giao thong, "Nhin lai hon 2 nam thu phi tu dong khong dung" [Toll lanes cut from 62 to 32 on the corridor after the ETC-only pilot], 2024. <https://www.baogiaothong.vn/nhin-lai-hon-2-nam-thu-phi-tu-dong-khong-dung-19224112407523824.htm> (accessed 10 October 2026) |
| Legal basis of electronic toll collection | Prime Minister of Vietnam, Decision 19/2020/QD-TTg of 17 June 2020 on non-stop electronic toll collection (ETC). <https://vanban.chinhphu.vn/default.aspx?pageid=27160&docid=200267> (accessed 10 October 2026) |

### Source documents as downloaded

The documents themselves are not redistributed here. Download them from the links above; their SHA-256
checksums (files as downloaded on 10 October 2026) are:

| Document | SHA-256 |
|---|---|
| Vietnam Road Administration spreadsheet (annex "Tong hop BOT 2025", sheet for 2024) | `e0110b9df4ca0138bb0b556c059821d254e990a9dca545d67612808d8f1570df` |
| VIDIFI audited consolidated financial statements 2024 (PDF) | `0ebe9c4af4a44db1307fb2dbbf9c46d3aa9cf545eb7023b628c3df87942b3480` |
| VIDIFI tariff notice, expressway fares by section (image) | `05a191950325f818bbe81f0c33239d786c393ee488f521774c04fb73193b241d` |
| VIDIFI tariff notice, QL5 fares per station (image) | `faca182938899fbef2f9a348dd54256d7a36db148dd107cde57aefca639d7ab8` |

## Episodes of toll avoidance and protest (`data/episodes.csv`)

Columns: `episode`; `partial_cut_pct`, the largest price concession offered before the episode ended (%);
`resolution`, how the episode ended (`admin_ban`, `shutdown`, `full_exempt_local`, `station_removal_buyout`,
`relocation`, `enforcement_abandoned`). Details are in SI, Section 9 and Table S10.

**DT391/QL5 (Hai Duong)**

- Ministry of Finance of Vietnam, Circular 153/2015/TT-BTC of 2 October 2015 on road-use fee rates at the two toll stations on National Road 5. <https://congbao.chinhphu.vn/van-ban/thong-tu-so-153-2015-tt-btc-17978.htm> (accessed 10 October 2026)
- VietnamPlus (Vietnam News Agency), "Phi Quoc lo 5 tang, tinh lo 391 oan minh ganh xe tai trong lon" [Hai Duong Department of Transport counts on DT391: containers 9 per day in 2015, 628 in the first quarter of 2016, 1,458 on 24 to 26 March 2016], 14 April 2016. <https://www.vietnamplus.vn/phi-quoc-lo-5-tang-tinh-lo-391-oan-minh-ganh-xe-tai-trong-lon-post381177.vnp> (accessed 10 October 2026)
- Thuong hieu va Cong luan, "Hai Duong: Tinh lo 391 keu cuu!" [Container traffic routed via QL10 during the resurfacing of QL5 in Hai Phong], 15 April 2016. <https://thuonghieucongluan.com.vn/hai-duong-tinh-lo-391-keu-cuu-a27373.html> (accessed 10 October 2026)
- Ministry of Finance of Vietnam, Circular 254/2016/TT-BTC amending Circular 153/2015/TT-BTC on the road-use fee at the two toll stations on National Road 5, effective 20 November 2016. <https://congbao.chinhphu.vn/van-ban/thong-tu-so-254-2016-tt-btc-21457.htm> (accessed 10 October 2026)
- BNews/TTXVN, "Hai Duong cam xe tai 4 truc tro len luu thong tren duong tinh 391" [45,725 four-axle vehicles/month in Dec 2016; total ban from 15 Jan 2017], 2017. <https://bnews.vn/hai-duong-cam-xe-tai-4-truc-tro-len-luu-thong-tren-duong-tinh-391/33121.html> (accessed 10 October 2026)

**Cai Lay (Tien Giang)**

- VnExpress, "Bo Giao thong quyet dinh giam phi BOT Cai Lay" [cuts of 22 to 33% by vehicle class from 21 Aug 2017, 35,000 to 25,000 VND for cars; full exemption for classes 1 and 2 of residents of four nearby communes], 16 August 2017. <https://vnexpress.net/bo-giao-thong-quyet-dinh-giam-phi-bot-cai-lay-3628222.html> (accessed 10 October 2026)
- Tuoi Tre, "Doanh nghiep dau tu du an BOT Cai Lay dung truoc nguy co pha san" [tolling from 1 Aug 2017 met with payment in small notes; suspended from 14 Dec 2017 for five years; accrued interest >500 bn VND], 5 October 2022. <https://tuoitre.vn/doanh-nghiep-dau-tu-du-an-bot-cai-lay-dung-truoc-nguy-co-pha-san-20221005144826807.htm> (accessed 10 October 2026)

**Ben Thuy 1 (Nghe An)**

- Infonet/Vietnamnet, "Ve BOT cau Ben Thuy giam 50%, dan van om boc tien le ra mua ve de phan doi" [50% cut for vehicles of the nearby districts from 15 Apr 2017; protests continued], 6 April 2017. <https://infonet.vietnamnet.vn/ve-bot-cau-ben-thuy-giam-50-dan-van-om-boc-tien-le-ra-mua-ve-de-phan-doi-247373.html> (accessed 10 October 2026)
- Nghe An Provincial Police portal, "Lam the mien phi qua cau Ben Thuy 1 nhu the nao?!" [the 50% cut dropped on 11 Apr 2017; vehicles of classes 1 and 2 of residents of Hong Linh, Nghi Xuan, Hung Nguyen and Vinh exempt from 24 Apr 2017], April 2017. <https://congan.nghean.gov.vn/cai-cach-hanh-chinh/201704/lam-the-mien-phi-qua-cau-ben-thuy-1-nhu-the-nao-732607/index.htm> (accessed 10 October 2026)

**T2/QL91 (Can Tho)**

- Tuoi Tre, "Chon phuong an bo thu phi BOT tram T2 quoc lo 91" [Ministry and investor agree to end collection at T2 and to buy back the QL91B section with budget funds], 11 November 2019. <https://tuoitre.vn/chon-phuong-an-bo-thu-phi-bot-tram-t2-quoc-lo-91-20191111174327829.htm> (accessed 10 October 2026)
- Tuoi Tre, "Xoa tram thu phi T1 va T2 tren quoc lo 91: can 1.879 ti dong tu ngan sach" [protests from Jul 2017; collection at T2 halted from 25 May 2019 after the Vam Cong bridge opened; project revenue at 50% (2020) and 36% (2021) of the financial plan; state buyout proposal], 1 July 2022. <https://tuoitre.vn/xoa-tram-thu-phi-t1-va-t2-tren-quoc-lo-91-can-1-879-ti-dong-tu-ngan-sach-20220701173225341.htm> (accessed 10 October 2026)

**Cau Rac (Ha Tinh)**

- VnExpress, "Tram thu phi Cau Rac dung ban ve" [protests from Apr 2016; full exemption for residents of Cam Xuyen and Ky Anh from 27 Apr 2016; collection suspended from 21 Feb 2019], 21 February 2019. <https://vnexpress.net/tram-thu-phi-cau-rac-dung-ban-ve-3884171.html> (accessed 10 October 2026)
- Tien Phong, "Vi sao tram BOT Cau Rac dot ngot dung thu phi?" [collection suspended from 21 Feb 2019 because the investment had largely been recovered; protests in 2016 to 2018], 20 February 2019. <https://tienphong.vn/vi-sao-tram-bot-cau-rac-dot-ngot-dung-thu-phi-post1092149.tpo> (accessed 10 October 2026)

**Tan De (Thai Binh)**

- Nhan Dan, "Chinh thuc thao do, di doi Tram BOT Tan De" [station dismantled 24 Jan 2019 and relocated to the Dong Hung bypass], 2019. <https://nhandan.vn/chinh-thuc-thao-do-di-doi-tram-bot-tan-de-post348068.html> (accessed 10 October 2026)

**Hac Tri (Phu Tho)**

- VnExpress, "Pho thu tuong yeu cau tinh lai phi BOT o cau Hac Tri" [the Deputy Prime Minister asks for the toll and the routing to be reconsidered; affected households already exempt], 12 July 2016. <https://vnexpress.net/pho-thu-tuong-yeu-cau-tinh-lai-phi-bot-o-cau-hac-tri-3434675.html> (accessed 10 October 2026)
- Bao Dau thau, "He luy Du an BOT cau Viet Tri moi" [toll revenue of about 100 m VND/day at the new bridge; cuts of 80% (cars under 7 seats, Bach Hac ward) and 60% (Song Lo commune and Viet Tri city) not accepted], 17 March 2016. <https://baodauthau.vn/he-luy-du-an-bot-cau-viet-tri-moi-post11253.html> (accessed 10 October 2026)
- VietnamPlus (Vietnam News Agency), "Phu Tho: Nguoi dan quay tram thu phi BOT phan doi viec ep di cau moi" [residents surround the Hac Tri toll station; a 72% cut for residents meeting three conditions], 4 March 2016. <https://www.vietnamplus.vn/phu-tho-nguoi-dan-quay-tram-thu-phi-bot-phan-doi-viec-ep-di-cau-moi-post374296.vnp> (accessed 10 October 2026)
- Bao Dien tu Chinh phu, "Sua cau Viet Tri, cho phep o to duoi 7 cho luu thong" [residents broke the concrete barriers to the old bridge from 1 Aug 2016 and most vehicles returned to it; the old bridge to be repaired and opened to cars under 7 seats], 17 August 2016. <https://baochinhphu.vn/sua-cau-viet-tri-cho-phep-o-to-duoi-7-cho-luu-thong-102207643.htm> (accessed 10 October 2026)
