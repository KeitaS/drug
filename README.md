# drug

```python
import crModel
growth, (cobra, ribo) = crModel.crModel({"DB01082": dose})
print(growth)
```

---  
## 研究の流れ
<font size="5px">2016/06/28</font>  
単剤と二剤の結果が出力できるようにした。  
drugsの中身を変え、それに合わせてcreateModel, runの引数等を変更した。  
P_in, P_outは計算から求めるようにし、createModel内で計算結果をdrugsにいれることにした。  
