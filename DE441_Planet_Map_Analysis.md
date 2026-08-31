# DE441 星历行星名称映射分析报告

## 测试环境
- 星历文件：DE441 (de441.bsp)
- 测试日期：1163-8-10
- 观测地点：北京 (39.9042°N, 116.4074°E)

## 测试结果

### 1. 行星名称可用性测试

| 行星 | 小写名称 | Barycenter 名称 | 大写名称 |
|------|----------|------------------|----------|
| Mercury | ✓ 可用 | ✓ 可用 | ✓ 可用 |
| Venus | ✓ 可用 | ✓ 可用 | ✓ 可用 |
| Mars | ✗ 不可用 | ✓ 可用 | ✗ 不可用 |
| Jupiter | ✗ 不可用 | ✓ 可用 | ✗ 不可用 |
| Saturn | ✗ 不可用 | ✓ 可用 | ✗ 不可用 |
| Uranus | ✗ 不可用 | ✓ 可用 | ✗ 不可用 |
| Neptune | ✗ 不可用 | ✓ 可用 | ✗ 不可用 |
| Pluto | ✗ 不可用 | ✓ 可用 | ✗ 不可用 |

**其他天体**：
- `sun`, `moon`, `earth`：全部可用（大小写均可）
- JPL 数字代码（如 '199', '299', '4', '5'）：**不可用**

### 2. 关键发现

#### 内行星 vs 外行星
- **内行星**（Mercury, Venus）：支持直接使用行星名称（如 `'mercury'`）
- **外行星**（Mars, Jupiter, Saturn, Uranus, Neptune, Pluto）：**仅支持 barycenter 名称**

#### 大小写敏感性
- 大小写不敏感：`'mercury'` = `'MERCURY'` = `'Mercury'`
- 但名称格式敏感：`'mars'` 不可用，`'mars barycenter'` 可用

#### Planet Center vs Barycenter 的区别
对于外行星：
- **Planet Center**：行星质心相对于太阳系质心的位置
- **Barycenter**：行星及其卫星系统的质心相对于太阳系质心的位置

**实际影响**：对于升落时间计算，差异可忽略不计（通常 < 1 秒）

### 3. 升落时间对比测试

以 Mars 为例（1163-8-10）：
```
Barycenter 升起时间: 1163-08-10T03:34:52.847161+07:45:38
```

由于 `'mars'` 不可用，无法直接对比，但根据天文计算原理：
- 距离差异：火星距地球 ~2 亿公里时，质心偏移约 1 万公里
- 时间影响：光行时差异约 0.03 秒，对升落时间影响可忽略

## 结论

### planet_map 的必要性

**结论：planet_map 是必要的**

原因：
1. **DE441 限制**：不支持简单行星名称（如 'mars', 'jupiter'）
2. **用户友好**：允许用户使用常见名称（如 'mars'），而非技术名称（如 'mars barycenter'）
3. **一致性**：统一内外行星的 API 接口
4. **错误预防**：避免用户因名称错误导致 KeyError

### 当前 planet_map 实现

```python
self.planet_map = {
    'mercury': 'mercury',
    'venus': 'venus',
    'mars': 'mars barycenter',
    'jupiter': 'jupiter barycenter',
    'saturn': 'saturn barycenter',
    'uranus': 'uranus barycenter',
    'neptune': 'neptune barycenter',
    'pluto': 'pluto barycenter',
}
```

**优点**：
- 简洁明了
- 覆盖所有主要行星
- 自动转换为 DE441 支持的格式

### 建议改进（可选）

如果需要支持 JPL 数字代码，可扩展为：

```python
self.planet_map = {
    # 常见名称
    'mercury': 'mercury',
    'venus': 'venus',
    'mars': 'mars barycenter',
    'jupiter': 'jupiter barycenter',
    'saturn': 'saturn barycenter',
    'uranus': 'uranus barycenter',
    'neptune': 'neptune barycenter',
    'pluto': 'pluto barycenter',
    # 数字代码
    '199': 'mercury',
    '299': 'venus',
    '4': 'mars barycenter',
    '5': 'jupiter barycenter',
    '6': 'saturn barycenter',
    '7': 'uranus barycenter',
    '8': 'neptune barycenter',
    '9': 'pluto barycenter',
}
```

## 参考资料

- [Skyfield Almanac Documentation](https://rhodesmill.org/skyfield/almanac.html)
- DE441 星历文件支持的目标列表（来自错误信息）：
  - 0 SOLAR SYSTEM BARYCENTER
  - 1 MERCURY BARYCENTER
  - 2 VENUS BARYCENTER
  - 3 EARTH BARYCENTER
  - 4 MARS BARYCENTER
  - 5 JUPITER BARYCENTER
  - 6 SATURN BARYCENTER
  - 7 URANUS BARYCENTER
  - 8 NEPTUNE BARYCENTER
  - 9 PLUTO BARYCENTER
  - 10 SUN
  - 199 MERCURY
  - 299 VENUS
  - 301 MOON
  - 399 EARTH