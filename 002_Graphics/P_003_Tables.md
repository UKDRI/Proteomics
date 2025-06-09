# 📊 Tables in Markdown

---

## Basic Table

Create a simple table using pipes (`|`) and hyphens (`-`):

```markdown
| Header 1   | Header 2   | Header 3   |
|------------|------------|------------|
| Row 1 Col 1| Row 1 Col 2| Row 1 Col 3|
| Row 2 Col 1| Row 2 Col 2| Row 2 Col 3|
```

| Header 1   | Header 2   | Header 3   |
|------------|------------|------------|
| Row 1 Col 1| Row 1 Col 2| Row 1 Col 3|
| Row 2 Col 1| Row 2 Col 2| Row 2 Col 3|

---

## Table with Alignment

```markdown
| Left-Aligned | Center-Aligned | Right-Aligned |
|:-------------|:--------------:|--------------:|
| Row 1 Col 1  | Row 1 Col 2    | Row 1 Col 3   |
| Row 2 Col 1  | Row 2 Col 2    | Row 2 Col 3   |
```

| Left-Aligned | Center-Aligned | Right-Aligned |
|:-------------|:--------------:|--------------:|
| Row 1 Col 1  | Row 1 Col 2    | Row 1 Col 3   |
| Row 2 Col 1  | Row 2 Col 2    | Row 2 Col 3   |

---


## Table with Multi-Line Text

```markdown
| Feature      | Description                          |
|--------------|--------------------------------------|
| Example        | Example Example Example Example |
| Example     | Example Example Example Example|
| Example  |   Example Example Example Example   |
```

| Feature      | Description                          |
|--------------|--------------------------------------|
| Example        | Example Example Example Example |
| Example     | Example Example Example Example|
| Example  |   Example Example Example Example   |

---

## Complex Tables with HTML

```html
<table>
  <tr>
    <th colspan="2">Example Table</th>
  </tr>
  <tr>
    <td>example</td>
    <td>example example</td>
  </tr>
  <tr>
    <td>example</td>
    <td>example example</td>
  </tr>
  <tr>
    <td colspan="2" align="center">example example example example</td>
  </tr>
</table>
```

<table>
  <tr>
    <th colspan="2">Example Table</th>
  </tr>
  <tr>
    <td>example</td>
    <td>example example</td>
  </tr>
  <tr>
    <td>example</td>
    <td>example example</td>
  </tr>
  <tr>
    <td colspan="2" align="center">example example example example</td>
  </tr>
</table>

---

## 🙌 Acknowledgments

📚 Use this space to list resources

* 📌 [Reference 1](https://example.com)
* 📌 [Reference 2](https://example.com)
* 📌 [Reference 3](https://example.com)