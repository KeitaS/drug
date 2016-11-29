# drug

```python
import crModel
growth, (cobra, ribo) = crModel.crModel({"DB01082": dose})
print(growth)
```

---
## ディレクトリの説明
- <font size="5px">**images**</font>  
  - <font size="4px">**2016**</font>  
  - <font size="4px">**ribo4**</font>
    - **heatmap_IC30.png**  
      IC30 * 2を両端に取ったヒートマップ
    - **heatmap_neweval_1.png**  
      新しい評価基準(1)を設けて，ヒートマップを作成(傾き１のみ)．  
      両端 IC30．
    - **heatmap_neweval_2.png**  
      新しい評価基準(1)を設けて，ヒートマップを作成(傾き１のみ)．  
      両端 IC30 * 2．
    - **heatmap_neweval_3.png**  
      新しい評価基準(1)を設けて，ヒートマップを作成(傾き[1/4, 1/2, 1, 2, 4])
      両端 IC30 * 2
    - **heatmap_neweval_4.png**  
      新しい評価基準(2)を設けて，ヒートマップを作成(傾き[1/4, 1/2, 1, 2, 4])
      両端 IC30 * 2
    - **heatmap_neweval_5.pnt**  
      新しい評価基準(2)を設けて，ヒートマップを作成
      1:1のときの中点をもとに，傾きに合わせて両端を作成できるようにした．
    - **linechart_neweval.png**
      新しい評価基準(1)をもちいて，様々な切片を割り振ったときのLineChart.

- <font size="5px">**metabolic_model**</font>  
  FBAを用いて作成した，代謝経路阻害剤モデル．

- <font size="5px">**model**</font>

- <font size="5px">**past-script**</font>  
  過去のriboモデル．

## riboモデルの説明
- <font size="4px">**ribo5.py**</font>  
  antagonisticが起こりうる2つのパターンを入れたモデル．  

- <font size="4px">**ribo6.py**</font>  
  synergisticの機構を加えたモデル．  
  機能型リボソームに結合して解離するようにモデリング．

---
## 内容
nature2006の評価基準を用いると，Additiveが表現できない．→ 新たな評価基準を設定．  
- <font size="4px">新しい評価基準(1)</font>  
  heatmap_IC30.pngをもとに，切片を段階的にとり，傾きを決めてその範囲でLineChartを描いたときにどのような波形になるかで，Additive, Synergistic, antagonisticを判断する．  
  単調増加，単調減少ならAdditive，上に凸ならAntagonistic，下に凸ならSynergistic．  

- <font size="4px">新しい評価基準(2)</font>  
  新しい評価基準(1)をベースに，Antagonistic,Synergisticになったときに，線形グラフからどれだけ離れているかで度合いを数値化．  
