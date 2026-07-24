#!/bin/bash

# 設定：只要有指令執行失敗 (回傳值非 0) 就立即中斷腳本
set -e
# 設定：在管線 (pipeline) 指令中，只要其中一個失敗就視為整體失敗
set -o pipefail

# 檢查使用者是否輸入了 3 個參數
if [ "$#" -ne 3 ]; then
    echo ""
    echo "用法: $0 <輸入的 VCF 檔案> <輸出的 VCF.gz 檔案> <Reference FASTA 檔案>"
    echo ""
    echo "範例: $0 input.vcf output_norm.vcf.gz ~/PROJECT/database/Homo_sapiens_assembly38.fasta"
    echo ""
    exit 1
fi

INPUT_VCF="$1"
OUTPUT_VCF="$2"
REF_FASTA="$3"

# ==========================================
# 建立輸出目錄 (如果不存在)
# ==========================================
OUTPUT_DIR=$(dirname "$OUTPUT_VCF")
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "--> 輸出目錄不存在，正在建立目錄: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
fi

# ==========================================
# 定義 AWK 處理函數 (將 FORMAT 的 AF/VAF 和 DP 提取至 INFO 的 TAF 和 TDP)
# 使用函數封裝可以避免在 if/else 中重複寫兩次相同的 AWK 程式碼
# ==========================================
add_taf_tdp() {
    awk '
    BEGIN { FS="\t"; OFS="\t" }
    /^##INFO/ && !done { 
        print "##INFO=<ID=TAF,Number=A,Type=Float,Description=\"Tumor alternate allele fraction for PCGR\">"; 
        print "##INFO=<ID=TDP,Number=1,Type=Integer,Description=\"Tumor specific read depth for PCGR\">"; 
        done=1 
    }
    /^#/ { print $0; next }
    {
        # 找出 FORMAT 欄位中 AF (或 VAF) 與 DP 的對應位置
        split($9, fmt, ":"); af_idx = 0; dp_idx = 0;
        for(i=1; i<=length(fmt); i++) { 
            # 這裡同時支援 AF 或 VAF 的抓取
            if(fmt[i] == "AF" || fmt[i] == "VAF") af_idx = i; 
            if(fmt[i] == "DP") dp_idx = i; 
        }
        
        # 抓出樣本對應的數值
        split($10, smp, ":"); 
        af_val = (af_idx > 0) ? smp[af_idx] : ""; 
        dp_val = (dp_idx > 0) ? smp[dp_idx] : "";
        
        # 組裝要插入 INFO 的字串
        info_add = "";
        if(af_val != "" && af_val != ".") { info_add = "TAF=" af_val ";"; }
        if(dp_val != "" && dp_val != ".") { info_add = info_add "TDP=" dp_val ";"; }
        
        # 塞入 INFO (第 8 欄) 的最前方
        if(info_add != "") { $8 = info_add $8; }
        print $0;
    }'
}

# ==========================================
# 主流程開始
# ==========================================
echo "正在處理並提取 TAF/TDP: $INPUT_VCF ..."

# 先把 header 結果存進變數（避開 set -o pipefail + grep -q 造成的 SIGPIPE 問題）
HEADER_AF=$(bcftools view -h "$INPUT_VCF" | grep '##FORMAT=<ID=AF' || true)

if [ -n "$HEADER_AF" ]; then
    echo "--> [發現 FORMAT/AF] 執行標準篩選與 TAF/TDP 提取..."
    
    # 流程 A：直接使用 FORMAT/AF 篩選 -> 正規化 -> 排序輸出文字(-Ov) -> AWK 加入 TAF/TDP -> 壓縮輸出
    bcftools view -i 'FILTER="PASS" && FORMAT/DP>=10 && FORMAT/AF>=0.05' "$INPUT_VCF" | \
    bcftools norm -f "$REF_FASTA" -m -both | \
    bcftools sort -Ov | \
    add_taf_tdp | \
    bgzip -c > "$OUTPUT_VCF"

else
    echo "--> [未發現 FORMAT/AF] 動態計算 FORMAT/VAF、篩選並提取 TAF/TDP..."
    
    # 流程 B：補上 VAF -> 以 VAF 篩選 -> 正規化 -> 排序輸出文字(-Ov) -> AWK 加入 TAF/TDP -> 壓縮輸出
    bcftools plugin fill-tags "$INPUT_VCF" -- -t FORMAT/VAF | \
    bcftools view -i 'FILTER="PASS" && FORMAT/DP>=10 && FORMAT/VAF>=0.05' | \
    bcftools norm -f "$REF_FASTA" -m -both | \
    bcftools sort -Ov | \
    add_taf_tdp | \
    bgzip -c > "$OUTPUT_VCF"
fi

# ==========================================
# 建立索引
# ==========================================
echo "--> 建立 Tabix 索引中..."
tabix -p vcf "$OUTPUT_VCF"

echo "--> $OUTPUT_VCF 處理完成！"
echo "------------------------------------------------"